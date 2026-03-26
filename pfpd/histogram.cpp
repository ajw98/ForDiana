#include "histogram.h"

///////////////// HISTOGRAMMING CLASS ////////////////////

#define GRIDTOL_EQUALMAGNITUDE 1e-6 // Is |k_i| == |k_j| ?
template<typename FieldType, typename FloatType>
histogram<FieldType, FloatType>::histogram(FFTlayout &rcFFTlayout, bool kspace/*=true*/)
{
  // For preparing the histogram:
  // - Sort through vector magnitudes of k grid
  // - Then loop through SORTED list and count how many PWs map to each |k|
  //     in doing so, record in a 2D array: PWindx=Map[histogramgridindx][entry#]
  // - Also record each UNIQUE |k| using a push_back
  // - Later, reduce entire field onto head CPU and bin the data through averaging over star members.
  std::vector<RealType> grid_vectormagnitudes(0);
  std::vector<RealType> vec(0);
  // Need global set of vector magnitudes...can't combine from individual CPUs unless we do a global sort but make a local map...
  for(UInt m=0; m<rcFFTlayout.getNPWglobal(); m++)
  {
    if(kspace)
      vec = rcFFTlayout.kvec_of_indx(m,false);
    else
      vec = rcFFTlayout.rvec_of_indx_0centred(m,false);
    FloatType modvec=0.0;
    for(UInt i=0; i<vec.size(); i++)
      modvec += vec[i]*vec[i];
    grid_vectormagnitudes.push_back(std::sqrt(modvec));
  }
  // Sort by |r| or |k| using Quicksort algorithm
  std::vector<UInt> indx(rcFFTlayout.getNPWglobal()+1); // indx is 1 based, so need N+1 members and ignore 0 element
  indexx(rcFFTlayout.getNPWglobal(),&(grid_vectormagnitudes[0]),&(indx[0]));
  // Loop over |r| or |k| list and store:
  // Unique |r,k| as abcissae
  // Mapping for how each FFT index goes to a single |r,k| index
  if(indx[1]!=0)
    codeerror_abort("Unexpected grid ordering in histogram setup",__FILE__,__LINE__);
  // Initialize maps/grids
  std::vector<UInt> degeneratelist(1,0);
  FloatType previous=grid_vectormagnitudes[indx[1]];
  m_adAbcissae.clear();
  m_adAbcissae.push_back(previous);
  // Loop through sorted list and build unique abcissae and grid mappings
  FloatType tol=GRIDTOL_EQUALMAGNITUDE;
  for(UInt n=1; n<rcFFTlayout.getNPWglobal(); n++)
  {
    FloatType current = grid_vectormagnitudes[indx[n+1]];
    if(previous < current-tol || previous > current+tol)
    {
      // NEW |k| or |r| entry
      m_adAbcissae.push_back(current);
      m_anMapper.push_back(degeneratelist); // Push the old degenerate list to the map BEFORE starting the new list
      degeneratelist.clear();
    }
    degeneratelist.push_back(indx[n+1]); // Update the list of degenerate points that map to the present |vec|
    previous = current;
  }
  m_anMapper.push_back(degeneratelist);
  // Allocate histogram data storage
  m_adData.clear();
  m_adData.resize(m_adAbcissae.size(),FieldType(0.0));
}

template<typename FieldType, typename FloatType>
void histogram<FieldType, FloatType>::mapfromfield(Field<FieldType> &infield, UInt NumSamplesInAverage)
{
  // Make histogram of angular average - fetch data to local node, CPU storage
  std::vector<FieldType> data(infield.getFFTlayout().getNPWglobal());
  infield.copytoarray(data);
  // Generate the histogram for angular average
  for(UInt mk=0; mk < m_adAbcissae.size(); mk++)
  {
    m_adData[mk] = FieldType(0.0);
    UInt nk = m_anMapper[mk].size();
    for(UInt n=0; n<nk; n++)
      m_adData[mk] += data[m_anMapper[mk][n]];
    m_adData[mk] *= 1.0/FloatType(nk*NumSamplesInAverage);
  }
}

template<typename FieldType, typename FloatType>
void histogram<FieldType,FloatType>::writehistogram(std::string const &filename) const
{
  std::ofstream histstream;
  if(m_adData.size()>0)
  {
    histstream.open(filename.c_str());
    for(UInt ix=1; ix<m_adAbcissae.size(); ix++) // Note: omit k=0
        histstream << m_adAbcissae[ix] << " " << m_adData[ix].real() << "  " << m_adData[ix].imag() << "\n";
    histstream.close();
  }
}

// Template specialization stubs
// Class specializations
template class histogram<std::complex<double>,double>;
//template class histogram<double,double>;
#ifdef ENABLESP
template class histogram<std::complex<float>,float>;
//template class histogram<float,float>;
#endif
