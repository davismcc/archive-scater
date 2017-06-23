#include "scater.h"
#include "c++/H5Cpp.h"

/********************* A rechunking class ************************/

template<typename T, bool use_size>
class rechunker {
public: 
    rechunker(const std::string& input_file, const std::string& input_data, 
              const std::string& output_file, const std::string& output_data,
              double longdim, bool byrow) : 
        ihfile(H5std_string(input_file), H5F_ACC_RDONLY),
        ihdata(ihfile.openDataSet(H5std_string(input_data))),
        ihspace(ihdata.getSpace()),
        HDT(ihdata.getDataType()),
        ohfile(H5std_string(output_file), H5F_ACC_RDWR)
    {
        // Setting up the input structures.
        if (ihspace.getSimpleExtentNdims()!=2){ 
            throw std::runtime_error("rechunking is not supported for arrays");
        }
        ihspace.getSimpleExtentDims(dims);
        inspace.setExtentSimple(2, dims);
       
        H5::DSetCreatPropList cparms = ihdata.getCreatePlist();
        cparms.getChunk(2, chunk_dims);

        // Setting up the storage structure. 
        hsize_t new_chunk_ncols=chunk_ncols(), new_chunk_nrows=chunk_nrows();
        if (byrow) {
            N=int(longdim/chunk_ncols() + 0.5);
            new_chunk_ncols=longdim;
            if (new_chunk_ncols > ncols()) { new_chunk_ncols=ncols(); }
        } else {
            N=int(longdim/chunk_nrows() + 0.5);
            new_chunk_nrows=longdim;
            if (new_chunk_nrows > nrows()) { new_chunk_nrows=nrows(); }
        }
        if (N <= 0) {
            throw std::runtime_error("number of chunks to load in per iteration should be positive");
        }
        storage_dims[0]=new_chunk_ncols;
        storage_dims[1]=new_chunk_nrows;
        storage_offset[0]=0;
        storage_offset[1]=0;
        store_space.setExtentSimple(2, storage_dims);
        const int datasize=(use_size ? HDT.getSize() : 1);
        cached.resize(storage_dims[0]*storage_dims[1]*datasize);
    
        // Setting up the output structures.       
        if (byrow) { 
            out_chunk_dims[0]=new_chunk_ncols;
            out_chunk_dims[1]=1;
        } else {
            out_chunk_dims[0]=1;
            out_chunk_dims[1]=new_chunk_nrows;
        }
        cparms.setChunk(2, out_chunk_dims);
        
        ohspace.setExtentSimple(2, dims);
        ohdata=ohfile.createDataSet(output_data, HDT, ohspace, cparms); 
        return;
    }

    void set_column_start() {     
        currentcol=0;
        out_ncols()=0;
        out_colpos()=0;
        storage_colpos()=0;
        counter=0;
        return;
    }

    bool advance_by_column (bool inner=true) {
        bool finished=false;
        in_colpos()=currentcol;
        if (currentcol + chunk_ncols() > ncols()) { 
            in_ncols()=ncols() - currentcol;
            finished=true;
        } else {
            in_ncols()=chunk_ncols();
        }
        if (inner) { 
            out_ncols() += in_ncols();
        } else {
            out_colpos() = currentcol;
            out_ncols() = in_ncols();
        }
        return finished;
    }
   
    void clean_up_column (bool inner=true) { 
        currentcol+=chunk_ncols();
        if (!inner) { return; }
        
        storage_colpos()+=chunk_ncols();
        ++counter;
        if (counter==N || currentcol>=ncols()) { 
            write_out_data();
            out_colpos()+=out_ncols();
            out_ncols()=0;
            counter=0;
        } 
        return;
    }

    void set_row_start() {     
        currentrow=0;
        out_nrows()=0;
        out_rowpos()=0;
        storage_rowpos()=0;
        counter=0;
        return;
    }

    bool advance_by_row (bool inner=true) {
        bool finished=false;
        in_rowpos()=currentrow;
        if (currentrow + chunk_nrows() > nrows()) { 
            in_nrows()=nrows() - currentrow;
            finished=true;
        } else {
            in_nrows()=chunk_nrows();
        }
        if (inner) {
            out_nrows() += in_nrows();
        } else {
            out_rowpos() = currentrow;
            out_nrows() = in_nrows();
        }
        return finished;
    }

    void clean_up_row (bool inner=true) { 
        currentrow+=chunk_nrows();
        if (!inner) { return; }
        
        storage_rowpos()+=chunk_nrows();
        ++counter;
        if (counter==N || currentrow>=nrows()) { 
            write_out_data();
            out_rowpos()+=out_nrows();
            out_nrows()=0;
            counter=0;
        } 
        return;
    }

    void read_in_data () {
        inspace.selectHyperslab(H5S_SELECT_SET, in_count, in_offset);
        store_space.selectHyperslab(H5S_SELECT_SET, in_count, storage_offset);
        ihdata.read(cached.data(), HDT, store_space, inspace);
        return;
    } 

    void write_out_data() {
        ohspace.selectHyperslab(H5S_SELECT_SET, out_count, out_offset);
        storage_colpos()=0;
        storage_rowpos()=0;
        store_space.selectHyperslab(H5S_SELECT_SET, out_count, storage_offset);
        ohdata.write(cached.data(), HDT, store_space, ohspace);
    }

private:
    H5::H5File ihfile;
    H5::DataSet ihdata;
    H5::DataSpace ihspace;
    const H5::DataType HDT;

    hsize_t dims[2];    
    hsize_t chunk_dims[2];
    
    H5::DataSpace inspace;
    hsize_t in_offset[2];
    hsize_t in_count[2];

    hsize_t storage_dims[2];
    hsize_t storage_offset[2];
    H5::DataSpace store_space;
    std::vector<T> cached;    

    H5::H5File ohfile;
    H5::DataSpace ohspace;
    H5::DataSet ohdata;
    
    hsize_t out_chunk_dims[2];
    hsize_t out_offset[2];
    hsize_t out_count[2];

    int N; // Store the number of chunks
    hsize_t currentcol, currentrow;
    int counter; 

    const hsize_t& ncols () const { return dims[0]; }
    const hsize_t& nrows () const { return dims[1]; }
    const hsize_t& chunk_ncols () { return chunk_dims[0]; }
    const hsize_t& chunk_nrows () { return chunk_dims[1]; }

    hsize_t& in_colpos () { return in_offset[0]; }
    hsize_t& in_rowpos () { return in_offset[1]; }
    hsize_t& in_ncols () { return in_count[0]; }
    hsize_t& in_nrows () { return in_count[1]; }

    const hsize_t& storage_ncol () const { return storage_dims[0]; }
    const hsize_t& storage_nrow () const { return storage_dims[1]; }
    hsize_t& storage_colpos () { return storage_offset[0]; }
    hsize_t& storage_rowpos () { return storage_offset[1]; } 

    hsize_t& out_chunk_ncols () { return out_chunk_dims[0]; }
    hsize_t& out_chunk_nrows () { return out_chunk_dims[1]; }

    hsize_t& out_colpos () { return out_offset[0]; }
    hsize_t& out_rowpos () { return out_offset[1]; }
    hsize_t& out_ncols () { return out_count[0]; }
    hsize_t& out_nrows () { return out_count[1]; }
};

/************************** Secondary templated functions *********************/

template <typename T, bool use_size> 
SEXP rechunk_by_row(Rcpp::StringVector ifile, Rcpp::StringVector idata, 
        Rcpp::StringVector ofile, Rcpp::StringVector odata, Rcpp::NumericVector nelements) {

    rechunker<T, use_size> repacker(Rcpp::as<std::string>(ifile[0]), Rcpp::as<std::string>(idata[0]),
            Rcpp::as<std::string>(ofile[0]), Rcpp::as<std::string>(odata[0]), 
            nelements[0], true);
    repacker.set_row_start();

    // Reading in entire chunks, and writing them back.
    while (1) { 
        bool finished=repacker.advance_by_row(false);
        repacker.set_column_start();

        while (1) { 
            bool inner_finished=repacker.advance_by_column(true);
            repacker.read_in_data();
            repacker.clean_up_column(true);
            if (inner_finished) { break; }        
        }

        repacker.clean_up_row(false);
        if (finished) { break; }
    }

    return R_NilValue;
}

template <typename T, bool use_size> 
SEXP rechunk_by_col(Rcpp::StringVector ifile, Rcpp::StringVector idata, 
        Rcpp::StringVector ofile, Rcpp::StringVector odata, Rcpp::NumericVector nelements) {

    rechunker<T, use_size> repacker(Rcpp::as<std::string>(ifile[0]), Rcpp::as<std::string>(idata[0]),
            Rcpp::as<std::string>(ofile[0]), Rcpp::as<std::string>(odata[0]), 
            nelements[0], false);
    repacker.set_column_start();

    // Reading in entire chunks, and writing them back.
    while (1) { 
        bool finished=repacker.advance_by_column(false);
        repacker.set_row_start();

        while (1) { 
            bool inner_finished=repacker.advance_by_row(true);
            repacker.read_in_data();
            repacker.clean_up_row(true);
            if (inner_finished) { break; }        
        }

        repacker.clean_up_column(false);
        if (finished) { break; }
    }

    return R_NilValue;
}

/************************** The actual R-visible functions *********************/

SEXP rechunk_matrix(SEXP inname, SEXP indata, SEXP intype, SEXP outname, SEXP outdata, SEXP longdim, SEXP byrow) {
    BEGIN_RCPP
    // Figuring out the type.
    Rcpp::StringVector type(intype);
    if (type.size()!=1) {
        throw std::runtime_error("type should be a string");
    }
    const std::string choice=Rcpp::as<std::string>(type[0]);

    // Figuring out the new chunk pattern.
    Rcpp::LogicalVector br(byrow);
    if (br.size()!=1) {
        throw std::runtime_error("byrow should be a logical scalar");
    }
    const bool BR=br[0];

    // Figuring out various other things.
    Rcpp::StringVector ifile(inname), idata(indata), ofile(outname), odata(outdata);
    if (ifile.size()!=1 || idata.size()!=1 || ofile.size()!=1 || odata.size()!=1) {
        throw std::runtime_error("file and dataset names must be strings");
    }
    Rcpp::NumericVector nelements(longdim);
    if (nelements.size()!=1) {
        throw std::runtime_error("chunk dimensions should be integer vectors of length 2");
    }

    // Dispatching.
    if (BR) {
        if (choice=="double") {
            return rechunk_by_row<double, false>(ifile, idata, outname, odata, nelements);
        } else if (choice=="integer" || choice=="logical") { 
            return rechunk_by_row<int, false>(ifile, idata, outname, odata, nelements);
        } else if (choice=="character") {
            return rechunk_by_row<char, true>(ifile, idata, outname, odata, nelements);
        }
    } else {
        if (choice=="double") {
            return rechunk_by_col<double, false>(ifile, idata, outname, odata, nelements);
        } else if (choice=="integer" || choice=="logical") { 
            return rechunk_by_col<int, false>(ifile, idata, outname, odata, nelements);
        } else if (choice=="character") {
            return rechunk_by_col<char, true>(ifile, idata, ofile, odata, nelements);
        }
    }
    throw std::runtime_error("unsupported data type");
    END_RCPP
}
