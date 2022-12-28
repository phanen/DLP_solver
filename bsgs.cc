#include <ctime>
#include <iostream>
#include <cmath>
#include <map>
#include <random>
using namespace std;

// Toy version of bsgs on elliptic curve(EC)
// Weierstrass form: 
//            y^2 = x^3 + ax + b mod p
// optional parameters: a, b, p


// alias for convinence...
using i32 = int;
using u32 = unsigned int;
using i64 = long long;
using u64 = unsigned long long; 

// type of finite field Zp  
using zp_t = u64; // considering extendable to NTL::ZZp
// type of point on toy curve
// always use point_t* (never point_t)
using point_t = pair<zp_t, zp_t>; 

ostream & operator<<(ostream & os, const point_t *P) {

  if (P)
    os << "(" << P->first << ", " << P->second << ")" << ends;
  else
    os << "(inf, inf)" << ends;
  return os;

}

// parameters of elliptic curve
struct param_t {
public:
  zp_t a, b, p; // enough for a toy

  param_t() {}
  param_t(zp_t a_, zp_t b_, zp_t p_) : a(a_), b(b_), p(p_)  {}

} ecc;  // a global public curve


// random number
mt19937 mt(time(nullptr));


// finite field (Zp) arithmetic
inline zp_t zp_add(zp_t l, zp_t r) { return (l + r) % ecc.p; }
inline zp_t zp_sub(zp_t l, zp_t r) { return (ecc.p + l - r) % ecc.p; }
inline zp_t zp_mul(zp_t l, zp_t r) { return (l * r) % ecc.p; }

inline zp_t zp_inv(zp_t a) {

  i64 t = 0, nt = 1;
  i64 r = ecc.p, nr = a;
  i64 tt, rr;

  while (nr) {
    auto mod = r / nr;
    tt = t; rr = r;
    t = nt; r = nr;
    nt = tt - mod * nt;
    nr = rr - mod * nr;
  }
  // if (r > 1) "not invertible";
  return t < 0 ? t + ecc.p : t;

}

bool ec_check(const point_t *P) {

  return
  zp_mul(P->second, P->second) == 
  zp_add( ecc.b, 
         zp_add(
           zp_mul(P->first, zp_mul(P->first, P->first)),
           zp_mul(ecc.a, P->first)
         )
         );

}

// R = P + P, assert P is not nullptr
inline void ec_dbl(const point_t *P, point_t *&R) { 
  
  if (!P) { delete R; R = nullptr; return; }
  
  auto xp = P->first, yp = P->second;
  if (!R) R = new point_t;

  auto lmd = zp_mul(
      zp_add(zp_mul(3, zp_mul(xp, xp)), ecc.a), 
      zp_inv(zp_mul(2, yp))
  );

  // cout << "lmd: " << lmd << endl;

  R->first = zp_sub(zp_mul(lmd, lmd), zp_add(xp, xp));
  R->second = zp_sub(zp_mul(lmd, zp_sub(xp, R->first)), yp);

}

// R = P + Q
inline void ec_add(const point_t *P, const point_t *Q, point_t *&R) { 

  // avoid alias reference
  point_t *PP = P ? new point_t(*P) : nullptr;
  point_t *QQ = Q ? new point_t(*Q) : nullptr;
  if (!R) R = new point_t;

  if (!PP) { 
    if (!QQ) { delete R; R = nullptr; return; }
    *R = *QQ;
    return;
  }
  if (!QQ) { 
    *R = *PP;
    return;
  }
  
  auto xp = PP->first, yp = PP->second;
  auto xq = QQ->first, yq = QQ->second;  
  
  if (xp == xq) {
    if (yp == yq) { ec_dbl(PP, R); }
    else { delete R; R = nullptr; }
    return;
  }
 
  auto lmd = zp_mul(zp_inv(zp_sub(xq, xp)), zp_sub(yq, yp));

  R->first = zp_sub(zp_mul(lmd, lmd), zp_add(xp, xq));
  R->second = zp_sub(zp_mul(lmd, zp_sub(xp, R->first)), yp);

}

// Q = -P
inline void ec_inv(const point_t *P, point_t *Q) {
  
  if (!P) { delete Q; Q = nullptr; return; }
  if (!Q) Q = new point_t;
  *Q = point_t(P->first, ecc.p - P->second); 

}

// Q = kP (double-and-add)
inline void ec_smul(const point_t *P, u64 k, point_t *&Q) {
  
  point_t *PP = new point_t(*P);
  u64 mask = 1ull << 31;
  if (Q) { delete Q; Q = nullptr; } 

  while (mask) {
    ec_dbl(Q, Q);
    if (k & mask) ec_add(Q, PP, Q);
    mask >>= 1;
  }

}


// return n (nP = Q)
// order ?? <= p + 1 + 2 * sqrt(p)
u64 bsgs(const point_t *P, const point_t *Q) {

  map<point_t, u64> tab;
  auto n = ecc.p + 1 + 2 * floor<u64>(sqrt(ecc.p));
  auto m = ceil<u64>(sqrt(n)) + 1;

  point_t *L = nullptr;
  for (u64 i = 1; i <= m; ++i) {
    ec_add(L, P, L); // baby step
    if (!L) continue;
    if (*L == *Q) return i;
    tab[*L] = i;
  } 

  point_t *R = nullptr, *G = nullptr;
  
  ec_smul(P, m, G);
  ec_inv(G, G);

  ec_add(Q, nullptr, R);

  for (u64 i = 1; i <= m; ++i) {

    ec_add(R, G, R); // gaint step
    if (!R) return i * m;
    if (tab.count(*R))  return tab[*R] + i * m;

  }

  return 0; // no solution
}

// return x^n
inline zp_t zp_pow(zp_t x, zp_t n) {

  u64 mask = 1ull << 31;
  zp_t y = 1; 
  while (mask) {
    y = zp_mul(y, y);
    if (n & mask) y = zp_mul(y, x);
    mask >>= 1;
  }
 
  return y;

}


// determine if x is a quadratic residue on Zp
inline bool zp_is_quad_res(zp_t x) {

  // auto mod8 = ecc.p % 8;
  // return (mod8 == 4 && x % 4 == 1 || mod8 == 0 && x % 8 == 1);
  return 1 == zp_pow(x, (ecc.p - 1) / 2);
}

// return y (y^2 = x) 
inline zp_t zp_clac_quad_res(zp_t x) {

  zp_t q = ecc.p - 1, s = 0;
  while (!(q & 1)) {
    q >>= 1; ++s;
  }  

  // search for a quadratic non-residue
  zp_t z;
  while (1) 
    if (!zp_is_quad_res(z = mt() % ecc.p)) break;  

  auto m = s;
  auto c = zp_pow(z, q);
  auto t = zp_pow(x, q);
  auto r = zp_pow(x, (q + 1) >> 1);

  while (1) {

    if (!t) return 0;
    if (t == 1) return r;
    u64 i = 0, ip = 1;

    while (i < m) {
      auto cur = zp_pow(t, ip);
      if (cur == 1) break;
      ++i;
      ip <<= 1;
    }

    auto b = zp_pow(c, 1 << (m - i - 1));
    m = i;
    c = zp_mul(b, b);
    t = zp_mul(t, c);
    r = zp_mul(r, b);
  }

  return 0;

}

void ec_gen_rand(point_t *P) {

  if (!P) P = new point_t;
  zp_t x, y, rhs;

  while (1) {
    
    x = mt() % ecc.p;
    rhs = zp_add(
        zp_add(zp_mul(x, zp_mul(x, x)), zp_mul(ecc.a, x)), ecc.b);
    // determine if rhs is a quadratic residue 
    if (zp_is_quad_res(rhs)) break;
  }
    // calculate y using tonelli-shanks algorithm
    // y^2 = rhs mod p
    y = zp_clac_quad_res(rhs);
    *P = {x, y};

}


u64 bf(const point_t *P, const point_t *Q) {

  point_t *T = nullptr;
  for (u64 i = 1; i < 20000; ++i) {
    ec_add(T, P, T);
    // cout << T << endl;
    if (!T) cout << i << endl;
    if (T && *T == *Q) return i;
  }
  return 0;

}

void bfm(const point_t *P, u64 k, point_t *&Q) {

  point_t *PP = new point_t(*P);
  if (Q) { delete Q; Q = nullptr; }
  for (u64 i = 1; i <= k; ++i) {
    ec_add(Q, PP, Q);
  }

}


void test_ec () {

  point_t *P = new point_t, *Q = new point_t, *T = new point_t;

  while (1) {

    ec_gen_rand(P);
    ec_gen_rand(Q);
 
    zp_t res = bsgs(P, Q);

    if (res == 0) continue; // no solution

    cout << "P: "<< P << endl;
    cout << "Q: "<< Q << endl;
    cout << "DLP(P, Q) = " << res << endl;

  }

}

zp_t zp_pow_bf (zp_t x, u64 n) {

  auto y = 1;
  while (n--)
    y = zp_mul(x, y); 
  return y;

}


// test zp
void test_zp() {

  zp_t a = 12;

  // test zp_pow
  for (u64 k = 1; k < 12345; ++k)
    if ((zp_pow_bf(a, k) != zp_pow(a, k))) cout << k << endl;

  // test zp_inv
  zp_t t = mt();

  for (u64 k = 1; k < 1234; ++k)
    if (t * zp_inv(t) % ecc.p != 1) cout << t << endl;

}


void ec_init() {

  ecc.a = 1253;
  ecc.b = 3271;
  ecc.p = 179321;

}

int main () {

  ec_init();
  // test_zp();
  test_ec();

}
