/***************************************************************
*
* Class defining a wrapper around the field class for better
* use of scratch space
* 
* DRT -- Wed, 19 Aug 2015
*
****************************************************************/

#ifndef _SMARTFIELD_H
#define _SMARTFIELD_H

#include "global.h"
#include "FieldStack.h"
#include "Field.h"

void init_SmartField_statics( UInt fd_size );
void delete_SmartField_statics();
void SmartField_unit_tests();
void test_time();

class SmartField : public FieldStack
{

  public:

    // --- constructor(s) and destructor ---

    SmartField();
    ~SmartField();

    SmartField( const SmartField &rhs ); // copy constructor w/ SmartField
    SmartField( const std::vector< Field<FieldType> > &rhs ); // copy constructor w/ vector of fields
    SmartField( const std::vector< std::vector<FieldType> > &rhs ); // copy constructor w/ vector of vectors
    SmartField( const Field<FieldType> &rhs ); // copy constructor w/ Field
    SmartField( const std::vector<FieldType> &rhs ); // copy constructor w/ STL vector
    SmartField( const FieldType rhs ); // copy constructor w/ constant

    void CreateField();

    // --- operators (field algebra) ---

    Field<FieldType> & operator[](const UInt i) { return (*data_ptr[i]); };
    const Field<FieldType> & operator[](const UInt i) const { return (*data_ptr[i]); };
    
    SmartField & operator=(const SmartField &rhs); // SmartField = SmartField
    SmartField & operator=(const std::vector< Field<FieldType> > &rhs); // SmartField = Vector of Fields
    SmartField & operator=(const std::vector< Field<RealType> > &rhs); // SmartField = Vector of Fields
    SmartField & operator=(const Field<FieldType> &rhs); // SmartField = Field
    SmartField & operator=(const std::vector< std::vector<FieldType> > &rhs); // SmartField = Vector of Vectors 
    SmartField & operator=(const std::vector<FieldType> &rhs); // Field = STL vector
    SmartField & operator=(const FieldType rhs); // Field = constant

    // --- Compound assignment ops ---
    // (Note: Kris' Field class provides assignments to,
    // but not compound operators with STL vectors.)

    SmartField & operator+=(const SmartField &rhs);
    SmartField & operator+=(const std::vector< Field<FieldType> > &rhs);
    SmartField & operator+=(const Field<FieldType> &rhs);
    SmartField & operator+=(FieldType rhs);

    SmartField & operator-=(const SmartField &rhs);
    SmartField & operator-=(const std::vector< Field<FieldType> > &rhs);
    SmartField & operator-=(const Field<FieldType> &rhs);
    SmartField & operator-=(const FieldType rhs);

    SmartField & operator*=(const SmartField &rhs);
    SmartField & operator*=(const std::vector< Field<FieldType> > &rhs);
    SmartField & operator*=(const Field<FieldType> &rhs);
    SmartField & operator*=(FieldType rhs);
    SmartField & operator*=(RealType rhs);

    SmartField & prodconjg(const SmartField &rhs);
    SmartField & prodconjg(const Field<FieldType> &rhs);

    SmartField & operator/=(const SmartField &rhs);
    SmartField & operator/=(const std::vector< Field<FieldType> > &rhs);
    SmartField & operator/=(const Field<FieldType> &rhs);
    SmartField & operator/=(const FieldType rhs);

    // --- member functions ---

    // like a subscript operator, [], but evaluates a smartfield at a vector of indices
    void subscript(const SmartField & rhs, const std::vector<UInt> idx);
    void slice(const SmartField & rhs, const std::vector<UInt> idx);

    // (mapping to Field class member functions)

    // initialization
    void fill(enum_fieldinitmethod filltype, Vec1dReal parameters=Vec1dReal(0), 
              bool screenoutputtype=false);
    void fillnoise(UInt noisetype);

    // fft operations
    void fft_rtok( bool applyscale=true );
    void fft_rtonegk( bool applyscale=true );
    void fft_ktor();
    inline void fft( bool applyscale=true ) { fft_rtok(applyscale);  };
    inline void ifft() { fft_ktor(); };
  
    // field manipulation/arithmetic setters
    void setaverage(const FieldType &value);
    void conjg();
    void zero();
    void zeroreal();
    void zeroimag();
    void zeropos();
    void zeroneg();
    void zero_k0(); // zero out the k=0 mode in the PS fields

    // data manipulation setters
    void setName(std::string szName_in);
    void setflag_inrealspace(bool bInRealSpace);
    void axpby_inplace(const SmartField &y, RealType a, RealType b);
    void axpby_inplace(const SmartField &y, Vec1dReal a, Vec1dReal b);
    void xpby_inplace(const SmartField &y, RealType b);
    void xpby_inplace(const SmartField &y, Vec1dReal b);
    void axpy_inplace(const SmartField &y, RealType a);
    void axpy_inplace(const SmartField &y, Vec1dReal a);
    //
    void axpby_inplace(const SmartField &y, FieldType a, FieldType b);
    void axpby_inplace(const SmartField &y, Vec1dFieldType a, Vec1dFieldType b);
    void xpby_inplace(const SmartField &y, FieldType b);
    void xpby_inplace(const SmartField &y, Vec1dFieldType b);
    void axpy_inplace(const SmartField &y, FieldType a);
    void axpy_inplace(const SmartField &y, Vec1dFieldType a);
    //
    void accumulateproduct_inplace(SmartField const &in1, SmartField const &in2);
    void accumulateproduct_inplace(SmartField const &in1, SmartField const &in2, double scale);
    //
    void exponentiate(const SmartField &in, RealType scale);
    void exponentiate(const SmartField &in);
    void sqrt(const SmartField &in);
    void log(const SmartField &in, RealType scale);
    void log(const SmartField &in);
    void logreal(const SmartField &in, RealType scale);
    void logreal(const SmartField &in);
    void setelement(const FieldType &value, UInt indx);

    // data manipulation getters
    bool getflag_inrealspace() const;
    const FFTlayout& getFFTlayout() const;
    FieldType getelement(UInt indx) const;
    void copytoarray(std::vector< std::vector<FieldType> > &out) const;
    void globalcopytoarray(std::vector< std::vector<FieldType> > &out) const;
    RealType l1norm() const;
    RealType l2norm() const;
    FieldType max() const;
    FieldType maxsigned() const;
    FieldType minsigned() const;
    UInt maxidx() const;
    FieldType getaverage() const;
    FieldType integrate() const;
    FieldType integrate_square() const;
    std::string getName() const;
    FieldType sumelem() const;
    FieldType sumelem_weighted(const SmartField &weights) const;
    FieldType sumelem_weighted(const Field<FieldType> &weights) const;

    // IO
    void writefield(std::string filename, bool omitimaginary);
    inline void writerealfield(std::string filename) { writefield(filename, true); };
    inline void writecplxfield(std::string filename) { writefield(filename, false); };
    void readfield(std::string filename);
    void debug_print_nelements(UInt num) const;

    // (other member functions)

    // directly obtain the field pointer to field "m"
    Field<FieldType>& GetField(UInt m = 0) const { return *(data_ptr[m]); };
    int GetSize() const { return data_ptr.size(); }

    // --- static member variables ---
    // * must assign before constructor is called!!
    static UInt BinaryInFlag;
    static UInt BinaryOutFlag;
    static UInt ICFlag;
    static UInt FDSize;
    static UInt NPW;

  private:

    // pointer to the vector of fields I'm creating
    std::vector< Field<FieldType>* > data_ptr; 

};

#endif //_SMARTFIELD_H
