

## ECDLP solver
A toy implementation of BSGS algorithm on elliptic curve

## Build
```
c++ bsgs.cc -o bsgs && ./bsgs
```

## Feature
- Support weierstrass form: $y^2 = x^3 + ax + b \bmod{p}$
- Optional parameter of ellitpic curve: a, b, p
- Automatically generate random point on curve (Legendre + T
onelli-Shanks algorithm)

## TODO
- Replace naive API with NTL API to support big integer.
- Other algrithm

## API
### Elliptic curve
```cpp
// parameters of elliptic curve
struct param_t {
public:
  zp_t a, b, p; // enough for a toy

  param_t() {}
  param_t(zp_t a_, zp_t b_, zp_t p_) : a(a_), b(b_), p(p_)  {}

} ecc;  // a global public curve
```

### Finite field (Zp) arithmetic
```cpp
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
```

### Elliptic curve arithmetic
```cpp
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

```

### BSGS
```cpp
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

```

## Reference

- https://en.wikipedia.org/wiki/Baby-step_giant-step
- https://crypto.stackexchange.com/questions/43614/how-to-find-the-generator-of-an-elliptic-curve
- https://en.wikipedia.org/wiki/Baby-step_giant-step
- https://en.wikipedia.org/wiki/Counting_points_on_elliptic_curves#Approaches_to_counting_points_on_elliptic_curves
- https://hyperelliptic.org/EFD/
