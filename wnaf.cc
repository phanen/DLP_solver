#include <iostream>
#include <pthread.h>
#include <vector>
using namespace std;


int to_digit(vector<int> naf) {
  
  int ret = 0;
  int base = 1;

  for (auto &b: naf) {
    ret += base * b;
    base <<= 1;
  }

  return ret;

}


int smod(int d, int mod) {
  int modd = mod >> 1;
  d %= mod;
  return d > modd ? d - mod : d; 
}
vector<int> to_naf(int d) {

  vector<int> ret;

  int w = 4;
  int mod = 1 << w; 

  int i = 0;
  while (d > 0) {

    if (d & 1) {
      ret.emplace_back(smod(d, mod));
      d -= ret[i];
    }
    else ret.emplace_back(0);
    
    d >>= 1;
    ++i;
  } 

  // auto mod2 = mod >> 1;
  // for (auto &each: ret) 
  //   if (each >= mod2) each -= mod;
  
  return ret;
}

int main() {

  vector<int> vd;
  int d = 100;

  for (int d = 0; d < 1000; ++d) {  

    vd = to_naf(d); 

    cout << d << endl;
    for (int i = vd.size() - 1; i >= 0; --i) 
      cout << vd[i] << ' ';
    cout << endl;
    // cout << to_digit(vd) << endl;

  }
}
