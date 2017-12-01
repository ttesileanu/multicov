#include "utils.h"

#include <stdexcept>

std::string to_string(const mxArray* v0)
{
  size_t len = mxGetM(v0)*mxGetN(v0);
  mxChar* v = mxGetChars(v0);

  // doing this manually to allow implicit conversions between mxChar and char
  std::string res(len, ' ');
  for (size_t i = 0; i < len; ++i)
    res[i] = v[i];

  return res;
}

Sequence to_seq(const mxChar* p, const std::vector<std::string>& alphabets,
    const std::vector<size_t> alphawidths)
{
  Sequence res;

  const size_t nalphas = alphabets.size();
  size_t crtidx = 0;
  for (size_t a = 0; a < nalphas; ++a) {
    const size_t w = alphawidths[a];
    const std::string& letters = alphabets[a];
    for (size_t i = 0; i < w; ++i) {
      const mxChar crtchar = p[crtidx];
      const size_t pos = letters.find(crtchar);
      if (pos == std::string::npos)
        throw std::runtime_error(std::string("Character ") + (char)crtchar +
          " not found in alphabet " + letters + ".");
      res.push_back(pos);

      ++crtidx;
    }
  }

  return res;
}
