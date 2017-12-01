#include "couplings.h"

#include <stdexcept>

void Couplings::set(
    const mxArray* dataorigin,
    const std::vector<std::string>& alphabets_,
    const std::vector<size_t>& alphawidths_)
{ 
  const size_t Jsize = mxGetM(dataorigin);
  data.setStorage(mxGetPr(dataorigin), Jsize, Jsize);

  size_t crtidx = 0;
  for (size_t a = 0; a < alphabets_.size(); ++a) {
    const std::string& letters = alphabets_[a];
    const size_t nletters_thisalpha = letters.length();
    const size_t w = alphawidths_[a];
    for (size_t i = 0; i < w; ++i) {
      nletters.push_back(nletters_thisalpha);
      startidx.push_back(crtidx);

      crtidx += nletters_thisalpha;
    }
  }

  if (crtidx != Jsize)
    throw std::runtime_error("Mismatch between coupling matrix and alphabet "
      "information.");
}
