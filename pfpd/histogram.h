#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include "Field.h"

template<typename FieldType, typename FloatType>
class histogram {
  public:
    std::vector<RealType> m_adAbcissae;
    std::vector<FieldType> m_adData;
    std::vector<std::vector<UInt> > m_anMapper; // Take FFT index, return ...

    // ctor
    histogram(FFTlayout &rcFFTlayout, bool kspace=true);
    // Methods
    void mapfromfield(Field<FieldType> &infield, UInt NumSamplesInAverage);
    // output
    void writehistogram(std::string const &filename) const;
};

#endif // _HISTOGRAM_H_
