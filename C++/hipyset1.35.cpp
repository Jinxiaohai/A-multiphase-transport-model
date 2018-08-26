#include <fem.hpp> // Fortran EMulation library of fable module

namespace AMPT {

using namespace fem::major_types;

void
acosh(...)
{
  throw std::runtime_error(
    "Missing function implementation: acosh");
}

void
asinh(...)
{
  throw std::runtime_error(
    "Missing function implementation: asinh");
}

void
dfour(...)
{
  throw std::runtime_error(
    "Missing function implementation: dfour");
}

void
digk(...)
{
  throw std::runtime_error(
    "Missing function implementation: digk");
}

void
djgk(...)
{
  throw std::runtime_error(
    "Missing function implementation: djgk");
}

std::complex<float>
fgk(...)
{
  throw std::runtime_error(
    "Missing function implementation: fgk");
}

void
four(...)
{
  throw std::runtime_error(
    "Missing function implementation: four");
}

void
hmeps(...)
{
  throw std::runtime_error(
    "Missing function implementation: hmeps");
}

void
ludbrb(...)
{
  throw std::runtime_error(
    "Missing function implementation: ludbrb");
}

void
pawt(...)
{
  throw std::runtime_error(
    "Missing function implementation: pawt");
}

struct common_ludat2
{
  arr<int, 2> kchg;
  arr<float, 2> pmas;
  arr<float> parf;
  arr<float, 2> vckm;

  common_ludat2() :
    kchg(dimension(500, 3), fem::fill0),
    pmas(dimension(500, 4), fem::fill0),
    parf(dimension(2000), fem::fill0),
    vckm(dimension(4, 4), fem::fill0)
  {}
};

struct common_ludat1
{
  arr<int> mstu;
  arr<float> paru;
  arr<int> mstj;
  arr<float> parj;

  common_ludat1() :
    mstu(dimension(200), fem::fill0),
    paru(dimension(200), fem::fill0),
    mstj(dimension(200), fem::fill0),
    parj(dimension(200), fem::fill0)
  {}
};

struct common_ludat4
{
  arr<fem::str<8> > chaf;

  common_ludat4() :
    chaf(dimension(500), fem::fill0)
  {}
};

struct common_lujets
{
  int n;
  arr<int, 2> k;
  arr<float, 2> p;
  arr<float, 2> v;

  common_lujets() :
    n(fem::int0),
    k(dimension(9000, 5), fem::fill0),
    p(dimension(9000, 5), fem::fill0),
    v(dimension(9000, 5), fem::fill0)
  {}
};

struct common_ludat3
{
  arr<int, 2> mdcy;
  arr<int, 2> mdme;
  arr<float> brat;
  arr<int, 2> kfdp;

  common_ludat3() :
    mdcy(dimension(500, 3), fem::fill0),
    mdme(dimension(2000, 2), fem::fill0),
    brat(dimension(2000), fem::fill0),
    kfdp(dimension(2000, 5), fem::fill0)
  {}
};

struct common_resdcy
{
  int nsav;
  int iksdcy;

  common_resdcy() :
    nsav(fem::int0),
    iksdcy(fem::int0)
  {}
};

struct common_pysubs
{
  int msel;
  arr<int> msub;
  arr<int, 2> kfin;
  arr<float> ckin;

  common_pysubs() :
    msel(fem::int0),
    msub(dimension(200), fem::fill0),
    kfin(dim1(2).dim2(-40, 40), fem::fill0),
    ckin(dimension(200), fem::fill0)
  {}
};

struct common_pypars
{
  arr<int> mstp;
  arr<float> parp;
  arr<int> msti;
  arr<float> pari;

  common_pypars() :
    mstp(dimension(200), fem::fill0),
    parp(dimension(200), fem::fill0),
    msti(dimension(200), fem::fill0),
    pari(dimension(200), fem::fill0)
  {}
};

struct common_pyint1
{
  arr<int> mint;
  arr<float> vint;

  common_pyint1() :
    mint(dimension(400), fem::fill0),
    vint(dimension(400), fem::fill0)
  {}
};

struct common_pyint4
{
  arr<float, 2> widp;
  arr<float, 2> wide;
  arr<float, 2> wids;

  common_pyint4() :
    widp(dim1(21, 40).dim2(0, 40), fem::fill0),
    wide(dim1(21, 40).dim2(0, 40), fem::fill0),
    wids(dim1(21, 40).dim2(3), fem::fill0)
  {}
};

struct common_pyint2
{
  arr<int> iset;
  arr<int, 2> kfpr;
  arr<float, 2> coef;
  arr<int, 3> icol;

  common_pyint2() :
    iset(dimension(200), fem::fill0),
    kfpr(dimension(200, 2), fem::fill0),
    coef(dimension(200, 20), fem::fill0),
    icol(dimension(40, 4, 2), fem::fill0)
  {}
};

struct common_pyint6
{
  arr<fem::str<28> > proc;

  common_pyint6() :
    proc(dim1(0, 200), fem::fill0)
  {}
};

struct common_pyint5
{
  arr<int, 2> ngen;
  arr<float, 2> xsec;

  common_pyint5() :
    ngen(dim1(0, 200).dim2(3), fem::fill0),
    xsec(dim1(0, 200).dim2(3), fem::fill0)
  {}
};

struct common_hparnt
{
  arr<float> hipr1;
  arr<int> ihpr2;
  arr<float> hint1;
  arr<int> ihnt2;

  common_hparnt() :
    hipr1(dimension(100), fem::fill0),
    ihpr2(dimension(50), fem::fill0),
    hint1(dimension(100), fem::fill0),
    ihnt2(dimension(50), fem::fill0)
  {}
};

struct common_hjcrdn
{
  arr<float, 2> yp;
  arr<float, 2> yt;

  common_hjcrdn() :
    yp(dimension(3, 300), fem::fill0),
    yt(dimension(3, 300), fem::fill0)
  {}
};

struct common_pyint3
{
  arr<float, 2> xsfx;
  arr<int, 2> isig;
  arr<float> sigh;

  common_pyint3() :
    xsfx(dim1(2).dim2(-40, 40), fem::fill0),
    isig(dimension(1000, 3), fem::fill0),
    sigh(dimension(1000), fem::fill0)
  {}
};

struct common_hstrng
{
  arr<int, 2> nfp;
  arr<float, 2> pphi;
  arr<int, 2> nft;
  arr<float, 2> pthi;

  common_hstrng() :
    nfp(dimension(300, 15), fem::fill0),
    pphi(dimension(300, 15), fem::fill0),
    nft(dimension(300, 15), fem::fill0),
    pthi(dimension(300, 15), fem::fill0)
  {}
};

struct common :
  fem::common,
  common_ludat2,
  common_ludat1,
  common_ludat4,
  common_lujets,
  common_ludat3,
  common_resdcy,
  common_pysubs,
  common_pypars,
  common_pyint1,
  common_pyint4,
  common_pyint2,
  common_pyint6,
  common_pyint5,
  common_hparnt,
  common_hjcrdn,
  common_pyint3,
  common_hstrng
{
  fem::variant_core common_ludatr;
  fem::cmn_sve rlu_sve;
  fem::cmn_sve lulist_sve;
  fem::cmn_sve ludecy_sve;
  fem::cmn_sve luboei_sve;
  fem::cmn_sve lugive_sve;
  fem::cmn_sve blockdata_ludata_sve;
  fem::cmn_sve pyinki_sve;
  fem::cmn_sve pyxtot_sve;
  fem::cmn_sve pygamm_sve;
  fem::cmn_sve pystfe_sve;
  fem::cmn_sve pystfu_sve;
  fem::cmn_sve pyspen_sve;
  fem::cmn_sve pymaxi_sve;
  fem::cmn_sve pyovly_sve;
  fem::cmn_sve pymult_sve;
  fem::cmn_sve pyinit_sve;
  fem::cmn_sve blockdata_pydata_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

//C
//C*********************************************************************
//C
int
lucomp(
  common& cmn,
  int const& kf)
{
  int return_value = fem::int0;
  // COMMON ludat2
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  //C
  //C...Purpose: to compress the standard KF codes for use in mass and decay
  //C...arrays; also to check whether a given code actually is defined.
  //C
  //C...Subdivide KF code into constituent pieces.
  return_value = 0;
  int kfa = fem::iabs(kf);
  int kfla = fem::mod(kfa / 1000, 10);
  int kflb = fem::mod(kfa / 100, 10);
  int kflc = fem::mod(kfa / 10, 10);
  int kfls = fem::mod(kfa, 10);
  int kflr = fem::mod(kfa / 10000, 10);
  //C
  //C...Simple cases: direct translation or special codes.
  if (kfa == 0 || kfa >= 100000) {
  }
  else if (kfa <= 100) {
    return_value = kfa;
    if (kf < 0 && kchg(kfa, 3) == 0) {
      return_value = 0;
    }
  }
  else if (kfls == 0) {
    if (kf == 130) {
      return_value = 221;
    }
    if (kf == 310) {
      return_value = 222;
    }
    if (kfa == 210) {
      return_value = 281;
    }
    if (kfa == 2110) {
      return_value = 282;
    }
    if (kfa == 2210) {
      return_value = 283;
    }
    //C
    //C...Mesons.
  }
  else if (kfa - 10000 * kflr < 1000) {
    if (kflb == 0 || kflb == 9 || kflc == 0 || kflc == 9) {
    }
    else if (kflb < kflc) {
    }
    else if (kf < 0 && kflb == kflc) {
    }
    else if (kflb == kflc) {
      if (kflr == 0 && kfls == 1) {
        return_value = 110 + kflb;
      }
      else if (kflr == 0 && kfls == 3) {
        return_value = 130 + kflb;
      }
      else if (kflr == 1 && kfls == 3) {
        return_value = 150 + kflb;
      }
      else if (kflr == 1 && kfls == 1) {
        return_value = 170 + kflb;
      }
      else if (kflr == 2 && kfls == 3) {
        return_value = 190 + kflb;
      }
      else if (kflr == 0 && kfls == 5) {
        return_value = 210 + kflb;
      }
    }
    else if (kflb <= 5 && kflc <= 3) {
      if (kflr == 0 && kfls == 1) {
        return_value = 100 + ((kflb - 1) * (kflb - 2)) / 2 + kflc;
      }
      else if (kflr == 0 && kfls == 3) {
        return_value = 120 + ((kflb - 1) * (kflb - 2)) / 2 + kflc;
      }
      else if (kflr == 1 && kfls == 3) {
        return_value = 140 + ((kflb - 1) * (kflb - 2)) / 2 + kflc;
      }
      else if (kflr == 1 && kfls == 1) {
        return_value = 160 + ((kflb - 1) * (kflb - 2)) / 2 + kflc;
      }
      else if (kflr == 2 && kfls == 3) {
        return_value = 180 + ((kflb - 1) * (kflb - 2)) / 2 + kflc;
      }
      else if (kflr == 0 && kfls == 5) {
        return_value = 200 + ((kflb - 1) * (kflb - 2)) / 2 + kflc;
      }
    }
    else if ((kfls == 1 && kflr <= 1) || (kfls == 3 && kflr <= 2) || (
      kfls == 5 && kflr == 0)) {
      return_value = 80 + kflb;
    }
    //C
    //C...Diquarks.
  }
  else if ((kflr == 0 || kflr == 1) && kflc == 0) {
    if (kfls != 1 && kfls != 3) {
    }
    else if (kfla == 9 || kflb == 0 || kflb == 9) {
    }
    else if (kfla < kflb) {
    }
    else if (kfls == 1 && kfla == kflb) {
    }
    else {
      return_value = 90;
    }
    //C
    //C...Spin 1/2 baryons.
  }
  else if (kflr == 0 && kfls == 2) {
    if (kfla == 9 || kflb == 0 || kflb == 9 || kflc == 9) {
    }
    else if (kfla <= kflc || kfla < kflb) {
    }
    else if (kfla >= 6 || kflb >= 4 || kflc >= 4) {
      return_value = 80 + kfla;
    }
    else if (kflb < kflc) {
      return_value = 300 + ((kfla + 1) * kfla * (kfla - 1)) / 6 + (
        kflc * (kflc - 1)) / 2 + kflb;
    }
    else {
      return_value = 330 + ((kfla + 1) * kfla * (kfla - 1)) / 6 + (
        kflb * (kflb - 1)) / 2 + kflc;
    }
    //C
    //C...Spin 3/2 baryons.
  }
  else if (kflr == 0 && kfls == 4) {
    if (kfla == 9 || kflb == 0 || kflb == 9 || kflc == 9) {
    }
    else if (kfla < kflb || kflb < kflc) {
    }
    else if (kfla >= 6 || kflb >= 4) {
      return_value = 80 + kfla;
    }
    else {
      return_value = 360 + ((kfla + 1) * kfla * (kfla - 1)) / 6 + (
        kflb * (kflb - 1)) / 2 + kflc;
    }
  }
  //C
  return return_value;
}

struct rlu_save
{
  fem::variant_bindings ludatr_bindings;
};

//C
//C*********************************************************************
//C
float
rlu(
  common& cmn,
  int const& /* idum */)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(rlu);
  common_variant ludatr(cmn.common_ludatr, sve.ludatr_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<int> mrlu(dimension(6));
      mbr<int> mrlu1;
      mbr<int> mrlu2;
      mbr<int> mrlu3;
      mbr<int> mrlu4;
      mbr<int> mrlu5;
      mbr<int> mrlu6;
      mbr<float> rrlu(dimension(100));
      mbr<float> rrlu98;
      mbr<float> rrlu99;
      mbr<float> rrlu00;
      ludatr.allocate(),
        equivalence(mrlu, mrlu1, mrlu2, mrlu3, mrlu4, mrlu5, mrlu6)
          .align<2>()
           .with<1>(arr_index(1))
          .align<3>()
           .with<1>(arr_index(2))
          .align<4>()
           .with<1>(arr_index(3))
          .align<5>()
           .with<1>(arr_index(4))
          .align<6>()
           .with<1>(arr_index(5))
          .align<7>()
           .with<1>(arr_index(6)),
        equivalence(rrlu, rrlu98, rrlu99, rrlu00)
          .align<2>()
           .with<1>(arr_index(98))
          .align<3>()
           .with<1>(arr_index(99))
          .align<4>()
           .with<1>(arr_index(100))
      ;
    }
  }
  /* arr_ref<int> mrlu( */ ludatr.bind<int>() /* , dimension(6)) */ ;
  int& mrlu1 = ludatr.bind<int>();
  int& mrlu2 = ludatr.bind<int>();
  int& mrlu3 = ludatr.bind<int>();
  int& mrlu4 = ludatr.bind<int>();
  int& mrlu5 = ludatr.bind<int>();
  /* int& mrlu6 */ ludatr.bind<int>();
  arr_ref<float> rrlu(ludatr.bind<float>(), dimension(100));
  float& rrlu98 = ludatr.bind<float>();
  float& rrlu99 = ludatr.bind<float>();
  float& rrlu00 = ludatr.bind<float>();
  int ij = fem::int0;
  int kl = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  int k = fem::int0;
  int l = fem::int0;
  int ii = fem::int0;
  float s = fem::float0;
  float t = fem::float0;
  int jj = fem::int0;
  int m = fem::int0;
  float twom24 = fem::float0;
  int i24 = fem::int0;
  float runi = fem::float0;
  //C
  //C...Purpose: to generate random numbers uniformly distributed between
  //C...0 and 1, excluding the endpoints.
  //C
  //C...Initialize generation from given seed.
  if (mrlu2 == 0) {
    ij = fem::mod(mrlu1 / 30082, 31329);
    kl = fem::mod(mrlu1, 30082);
    i = fem::mod(ij / 177, 177) + 2;
    j = fem::mod(ij, 177) + 2;
    k = fem::mod(kl / 169, 178) + 1;
    l = fem::mod(kl, 169);
    FEM_DO_SAFE(ii, 1, 97) {
      s = 0.f;
      t = 0.5f;
      FEM_DO_SAFE(jj, 1, 24) {
        m = fem::mod(fem::mod(i * j, 179) * k, 179);
        i = j;
        j = k;
        k = m;
        l = fem::mod(53 * l + 1, 169);
        if (fem::mod(l * m, 64) >= 32) {
          s += t;
        }
        t = 0.5f * t;
      }
      rrlu(ii) = s;
    }
    twom24 = 1.f;
    FEM_DO_SAFE(i24, 1, 24) {
      twom24 = 0.5f * twom24;
    }
    rrlu98 = 362436.f * twom24;
    rrlu99 = 7654321.f * twom24;
    rrlu00 = 16777213.f * twom24;
    mrlu2 = 1;
    mrlu3 = 0;
    mrlu4 = 97;
    mrlu5 = 33;
  }
  //C
  //C...Generate next random number.
  statement_130:
  runi = rrlu(mrlu4) - rrlu(mrlu5);
  if (runi < 0.f) {
    runi += 1.f;
  }
  rrlu(mrlu4) = runi;
  mrlu4 = mrlu4 - 1;
  if (mrlu4 == 0) {
    mrlu4 = 97;
  }
  mrlu5 = mrlu5 - 1;
  if (mrlu5 == 0) {
    mrlu5 = 97;
  }
  rrlu98 = rrlu98 - rrlu99;
  if (rrlu98 < 0.f) {
    rrlu98 += rrlu00;
  }
  runi = runi - rrlu98;
  if (runi < 0.f) {
    runi += 1.f;
  }
  if (runi <= 0 || runi >= 1.f) {
    goto statement_130;
  }
  //C
  //C...Update counters. Random number to output.
  mrlu3++;
  if (mrlu3 == 1000000000) {
    mrlu2++;
    mrlu3 = 0;
  }
  return_value = runi;
  //C
  return return_value;
}

//C
//C*********************************************************************
//C
float
ulmass(
  common& cmn,
  int const& kf)
{
  float return_value = fem::float0;
  // COMMON ludat1
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  // COMMON ludat2
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_ref<float> parf(cmn.parf, dimension(2000));
  //
  //C
  //C...Purpose: to give the mass of a particle/parton.
  //C
  //C...Reset variables. Compressed code.
  return_value = 0.f;
  int kfa = fem::iabs(kf);
  int kc = lucomp(cmn, kf);
  if (kc == 0) {
    return return_value;
  }
  parf(106) = pmas(6, 1);
  parf(107) = pmas(7, 1);
  parf(108) = pmas(8, 1);
  //C
  //C...Guarantee use of constituent masses for internal checks.
  int kfla = fem::int0;
  int kflb = fem::int0;
  int kflc = fem::int0;
  int kfls = fem::int0;
  int kflr = fem::int0;
  float pma = fem::float0;
  float pmb = fem::float0;
  float pmc = fem::float0;
  float pmspl = fem::float0;
  int kmul = fem::int0;
  if ((mstj(93) == 1 || mstj(93) == 2) && kfa <= 10) {
    return_value = parf(100 + kfa);
    if (mstj(93) == 2) {
      return_value = fem::max(0.f, return_value - parf(121));
    }
    //C
    //C...Masses that can be read directly off table.
  }
  else if (kfa <= 100 || kc <= 80 || kc > 100) {
    return_value = pmas(kc, 1);
    //C
    //C...Find constituent partons and their masses.
  }
  else {
    kfla = fem::mod(kfa / 1000, 10);
    kflb = fem::mod(kfa / 100, 10);
    kflc = fem::mod(kfa / 10, 10);
    kfls = fem::mod(kfa, 10);
    kflr = fem::mod(kfa / 10000, 10);
    pma = parf(100 + kfla);
    pmb = parf(100 + kflb);
    pmc = parf(100 + kflc);
    //C
    //C...Construct masses for various meson, diquark and baryon cases.
    if (kfla == 0 && kflr == 0 && kfls <= 3) {
      if (kfls == 1) {
        pmspl = -3.f / (pmb * pmc);
      }
      if (kfls >= 3) {
        pmspl = 1.f / (pmb * pmc);
      }
      return_value = parf(111) + pmb + pmc + parf(113) * fem::pow2(
        parf(101)) * pmspl;
    }
    else if (kfla == 0) {
      kmul = 2;
      if (kfls == 1) {
        kmul = 3;
      }
      if (kflr == 2) {
        kmul = 4;
      }
      if (kfls == 5) {
        kmul = 5;
      }
      return_value = parf(113 + kmul) + pmb + pmc;
    }
    else if (kflc == 0) {
      if (kfls == 1) {
        pmspl = -3.f / (pma * pmb);
      }
      if (kfls == 3) {
        pmspl = 1.f / (pma * pmb);
      }
      return_value = 2.f * parf(112) / 3.f + pma + pmb + parf(114) *
        fem::pow2(parf(101)) * pmspl;
      if (mstj(93) == 1) {
        return_value = pma + pmb;
      }
      if (mstj(93) == 2) {
        return_value = fem::max(0.f, return_value - parf(122) - 2.f *
          parf(112) / 3.f);
      }
    }
    else {
      if (kfls == 2 && kfla == kflb) {
        pmspl = 1.f / (pma * pmb) - 2.f / (pma * pmc) - 2.f / (pmb * pmc);
      }
      else if (kfls == 2 && kflb >= kflc) {
        pmspl = -2.f / (pma * pmb) - 2.f / (pma * pmc) + 1.f / (pmb * pmc);
      }
      else if (kfls == 2) {
        pmspl = -3.f / (pmb * pmc);
      }
      else {
        pmspl = 1.f / (pma * pmb) + 1.f / (pma * pmc) + 1.f / (pmb * pmc);
      }
      return_value = parf(112) + pma + pmb + pmc + parf(114) *
        fem::pow2(parf(101)) * pmspl;
    }
  }
  //C
  //C...Optional mass broadening according to truncated Breit-Wigner
  //C...(either in m or in m^2).
  float pm0 = fem::float0;
  float pmlow = fem::float0;
  float pmupp = fem::float0;
  if (mstj(24) >= 1 && pmas(kc, 2) > 1e-4f) {
    if (mstj(24) == 1 || (mstj(24) == 2 && kfa > 100)) {
      return_value += 0.5f * pmas(kc, 2) * fem::tan((2.f * rlu(cmn,
        0) - 1.f) * fem::atan(2.f * pmas(kc, 3) / pmas(kc, 2)));
    }
    else {
      pm0 = return_value;
      pmlow = fem::atan((fem::pow2(fem::max(0.f, pm0 - pmas(kc,
        3))) - fem::pow2(pm0)) / (pm0 * pmas(kc, 2)));
      pmupp = fem::atan(fem::pow2((pm0 + pmas(kc, 3))) - fem::pow2(
        pm0)) / (pm0 * pmas(kc, 2));
      return_value = fem::sqrt(fem::max(0.f, fem::pow2(pm0) + pm0 * pmas(kc,
        2) * fem::tan(pmlow + (pmupp - pmlow) * rlu(cmn, 0))));
    }
  }
  mstj(93) = 0;
  //C
  return return_value;
}

//C
//C*********************************************************************
//C
int
luchge(
  common& cmn,
  int const& kf)
{
  int return_value = fem::int0;
  // COMMON ludat2
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  //C
  //C...Purpose: to give three times the charge for a particle/parton.
  //C
  //C...Initial values. Simple case of direct readout.
  return_value = 0;
  int kfa = fem::iabs(kf);
  int kc = lucomp(cmn, kfa);
  if (kc == 0) {
  }
  else if (kfa <= 100 || kc <= 80 || kc > 100) {
    return_value = kchg(kc, 1);
    //C
    //C...Construction from quark content for heavy meson, diquark, baryon.
  }
  else if (fem::mod(kfa / 1000, 10) == 0) {
    return_value = (kchg(fem::mod(kfa / 100, 10), 1) - kchg(fem::mod(kfa / 10,
      10), 1)) * fem::pow((-1), fem::mod(kfa / 100, 10));
  }
  else if (fem::mod(kfa / 10, 10) == 0) {
    return_value = kchg(fem::mod(kfa / 1000, 10), 1) + kchg(fem::mod(kfa / 100,
      10), 1);
  }
  else {
    return_value = kchg(fem::mod(kfa / 1000, 10), 1) + kchg(fem::mod(kfa / 100,
      10), 1) + kchg(fem::mod(kfa / 10, 10), 1);
  }
  //C
  //C...Add on correct sign.
  return_value = return_value * fem::isign(1, kf);
  //C
  return return_value;
}

//C
//C*********************************************************************
//C
void
luname(
  common& cmn,
  int const& kf,
  str_ref chau)
{
  // COMMON ludat1
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  // COMMON ludat4
  str_arr_cref<1> chaf(cmn.chaf, dimension(500));
  //
  //C
  //C...Purpose: to give the particle/parton name as a character string.
  //C
  //C...Initial values. Charge. Subdivide code.
  chau = " ";
  int kfa = fem::iabs(kf);
  int kc = lucomp(cmn, kf);
  if (kc == 0) {
    return;
  }
  int kq = luchge(cmn, kf);
  int kfla = fem::mod(kfa / 1000, 10);
  int kflb = fem::mod(kfa / 100, 10);
  int kflc = fem::mod(kfa / 10, 10);
  int kfls = fem::mod(kfa, 10);
  int kflr = fem::mod(kfa / 10000, 10);
  //C
  //C...Read out root name and spin for simple particle.
  int len = fem::int0;
  int lem = fem::int0;
  if (kfa <= 100 || (kfa > 100 && kc > 100)) {
    chau = chaf(kc);
    len = 0;
    FEM_DO_SAFE(lem, 1, 8) {
      if (chau(lem, lem) != " ") {
        len = lem;
      }
    }
    //C
    //C...Construct root name for diquark. Add on spin.
  }
  else if (kflc == 0) {
    chau(1, 2) = chaf(kfla)(1, 1) + chaf(kflb)(1, 1);
    if (kfls == 1) {
      chau(3, 4) = "_0";
    }
    if (kfls == 3) {
      chau(3, 4) = "_1";
    }
    len = 4;
    //C
    //C...Construct root name for heavy meson. Add on spin and heavy flavour.
  }
  else if (kfla == 0) {
    if (kflb == 5) {
      chau(1, 1) = "B";
    }
    if (kflb == 6) {
      chau(1, 1) = "T";
    }
    if (kflb == 7) {
      chau(1, 1) = "L";
    }
    if (kflb == 8) {
      chau(1, 1) = "H";
    }
    len = 1;
    if (kflr == 0 && kfls == 1) {
    }
    else if (kflr == 0 && kfls == 3) {
      chau(2, 2) = "*";
      len = 2;
    }
    else if (kflr == 1 && kfls == 3) {
      chau(2, 3) = "_1";
      len = 3;
    }
    else if (kflr == 1 && kfls == 1) {
      chau(2, 4) = "*_0";
      len = 4;
    }
    else if (kflr == 2) {
      chau(2, 4) = "*_1";
      len = 4;
    }
    else if (kfls == 5) {
      chau(2, 4) = "*_2";
      len = 4;
    }
    if (kflc >= 3 && kflr == 0 && kfls <= 3) {
      chau(len + 1, len + 2) = "_" + chaf(kflc)(1, 1);
      len += 2;
    }
    else if (kflc >= 3) {
      chau(len + 1, len + 1) = chaf(kflc)(1, 1);
      len++;
    }
    //C
    //C...Construct root name and spin for heavy baryon.
  }
  else {
    if (kflb <= 2 && kflc <= 2) {
      chau = "Sigma ";
      if (kflc > kflb) {
        chau = "Lambda";
      }
      if (kfls == 4) {
        chau = "Sigma*";
      }
      len = 5;
      if (chau(6, 6) != " ") {
        len = 6;
      }
    }
    else if (kflb <= 2 || kflc <= 2) {
      chau = "Xi ";
      if (kfla > kflb && kflb > kflc) {
        chau = "Xi'";
      }
      if (kfls == 4) {
        chau = "Xi*";
      }
      len = 2;
      if (chau(3, 3) != " ") {
        len = 3;
      }
    }
    else {
      chau = "Omega ";
      if (kfla > kflb && kflb > kflc) {
        chau = "Omega'";
      }
      if (kfls == 4) {
        chau = "Omega*";
      }
      len = 5;
      if (chau(6, 6) != " ") {
        len = 6;
      }
    }
    //C
    //C...Add on heavy flavour content for heavy baryon.
    chau(len + 1, len + 2) = "_" + chaf(kfla)(1, 1);
    len += 2;
    if (kflb >= kflc && kflc >= 4) {
      chau(len + 1, len + 2) = chaf(kflb)(1, 1) + chaf(kflc)(1, 1);
      len += 2;
    }
    else if (kflb >= kflc && kflb >= 4) {
      chau(len + 1, len + 1) = chaf(kflb)(1, 1);
      len++;
    }
    else if (kflc > kflb && kflb >= 4) {
      chau(len + 1, len + 2) = chaf(kflc)(1, 1) + chaf(kflb)(1, 1);
      len += 2;
    }
    else if (kflc > kflb && kflc >= 4) {
      chau(len + 1, len + 1) = chaf(kflc)(1, 1);
      len++;
    }
  }
  //C
  //C...Add on bar sign for antiparticle (where necessary).
  if (kf > 0 || len == 0) {
  }
  else if (kfa > 10 && kfa <= 40 && kq != 0) {
  }
  else if (kfa == 89 || (kfa >= 91 && kfa <= 99)) {
  }
  else if (kfa > 100 && kfla == 0 && kq != 0) {
  }
  else if (mstu(15) <= 1) {
    chau(len + 1, len + 1) = "~";
    len++;
  }
  else {
    chau(len + 1, len + 3) = "bar";
    len += 3;
  }
  //C
  //C...Add on charge where applicable (conventional cases skipped).
  if (kq == 6) {
    chau(len + 1, len + 2) = "++";
  }
  if (kq ==  - 6) {
    chau(len + 1, len + 2) = "--";
  }
  if (kq == 3) {
    chau(len + 1, len + 1) = "+";
  }
  if (kq ==  - 3) {
    chau(len + 1, len + 1) = "-";
  }
  if (kq == 0 && (kfa <= 22 || len == 0)) {
  }
  else if (kq == 0 && (kfa >= 81 && kfa <= 100)) {
  }
  else if (kfa > 100 && kfla == 0 && kflb == kflc && kflb != 1) {
  }
  else if (kq == 0) {
    chau(len + 1, len + 1) = "0";
  }
  //C
}

//C
//C*********************************************************************
//C
float
ulangl(
  common& cmn,
  float const& x,
  float const& y)
{
  float return_value = fem::float0;
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  //
  //C
  //C...Purpose: to reconstruct an angle from given x and y coordinates.
  //C
  return_value = 0.f;
  float r = fem::sqrt(fem::pow2(x) + fem::pow2(y));
  if (r < 1e-20f) {
    return return_value;
  }
  if (fem::abs(x) / r < 0.8f) {
    return_value = fem::sign(fem::acos(x / r), y);
  }
  else {
    return_value = fem::asin(y / r);
    if (x < 0.f && return_value >= 0.f) {
      return_value = paru(1) - return_value;
    }
    else if (x < 0.f) {
      return_value = -paru(1) - return_value;
    }
  }
  //C
  return return_value;
}

//C
//C*********************************************************************
//C
float
plu(
  common& cmn,
  int const& i,
  int const& j)
{
  float return_value = fem::float0;
  // COMMON lujets
  int& n = cmn.n;
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_cref<float, 2> p(cmn.p, dimension(9000, 5));
  // COMMON ludat1
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  //
  //C
  //C...Purpose: to provide various real-valued event related data.
  //C
  //C...Set default value. For I = 0 sum of momenta or charges,
  //C...or invariant mass of system.
  return_value = 0.f;
  int i1 = fem::int0;
  int j1 = fem::int0;
  arr_1d<4, float> psum(fem::fill0);
  float pmr = fem::float0;
  float pr = fem::float0;
  if (i < 0 || i > mstu(4) || j <= 0) {
  }
  else if (i == 0 && j <= 4) {
    FEM_DO_SAFE(i1, 1, n) {
      if (k(i1, 1) > 0 && k(i1, 1) <= 10) {
        return_value += p(i1, j);
      }
    }
  }
  else if (i == 0 && j == 5) {
    FEM_DO_SAFE(j1, 1, 4) {
      psum(j1) = 0.f;
      FEM_DO_SAFE(i1, 1, n) {
        if (k(i1, 1) > 0 && k(i1, 1) <= 10) {
          psum(j1) += p(i1, j1);
        }
      }
    }
    return_value = fem::sqrt(fem::max(0.f, fem::pow2(psum(4)) -
      fem::pow2(psum(1)) - fem::pow2(psum(2)) - fem::pow2(psum(3))));
  }
  else if (i == 0 && j == 6) {
    FEM_DO_SAFE(i1, 1, n) {
      if (k(i1, 1) > 0 && k(i1, 1) <= 10) {
        return_value += luchge(cmn, k(i1, 2)) / 3.f;
      }
    }
  }
  else if (i == 0) {
    //C
    //C...Direct readout of P matrix.
  }
  else if (j <= 5) {
    return_value = p(i, j);
    //C
    //C...Charge, total momentum, transverse momentum, transverse mass.
  }
  else if (j <= 12) {
    if (j == 6) {
      return_value = luchge(cmn, k(i, 2)) / 3.f;
    }
    if (j == 7 || j == 8) {
      return_value = fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) + fem::pow2(p(i,
        3));
    }
    if (j == 9 || j == 10) {
      return_value = fem::pow2(p(i, 1)) + fem::pow2(p(i, 2));
    }
    if (j == 11 || j == 12) {
      return_value = fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) + fem::pow2(p(i,
        2));
    }
    if (j == 8 || j == 10 || j == 12) {
      return_value = fem::sqrt(return_value);
    }
    //C
    //C...Theta and phi angle in radians or degrees.
  }
  else if (j <= 16) {
    if (j <= 14) {
      return_value = ulangl(cmn, p(i, 3), fem::sqrt(fem::pow2(p(i,
        1)) + fem::pow2(p(i, 2))));
    }
    if (j >= 15) {
      return_value = ulangl(cmn, p(i, 1), p(i, 2));
    }
    if (j == 14 || j == 16) {
      return_value = return_value * 180.f / paru(1);
    }
    //C
    //C...True rapidity, rapidity with pion mass, pseudorapidity.
  }
  else if (j <= 19) {
    pmr = 0.f;
    if (j == 17) {
      pmr = p(i, 5);
    }
    if (j == 18) {
      pmr = ulmass(cmn, 211);
    }
    pr = fem::max(1e-20f, fem::pow2(pmr) + fem::pow2(p(i, 1)) + fem::pow2(p(i,
      2)));
    return_value = fem::sign(fem::log(fem::min((fem::sqrt(pr + fem::pow2(p(i,
      3))) + fem::abs(p(i, 3))) / fem::sqrt(pr), 1e20f)), p(i, 3));
    //C
    //C...Energy and momentum fractions (only to be used in CM frame).
  }
  else if (j <= 25) {
    if (j == 20) {
      return_value = 2.f * fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i,
        2)) + fem::pow2(p(i, 3))) / paru(21);
    }
    if (j == 21) {
      return_value = 2.f * p(i, 3) / paru(21);
    }
    if (j == 22) {
      return_value = 2.f * fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i,
        2))) / paru(21);
    }
    if (j == 23) {
      return_value = 2.f * p(i, 4) / paru(21);
    }
    if (j == 24) {
      return_value = (p(i, 4) + p(i, 3)) / paru(21);
    }
    if (j == 25) {
      return_value = (p(i, 4) - p(i, 3)) / paru(21);
    }
  }
  //C
  return return_value;
}

struct lulist_save
{
  arr<fem::str<4> > chdl;
  arr<fem::str<3> > chmo;

  lulist_save() :
    chdl(dimension(7), fem::fill0),
    chmo(dimension(12), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
lulist(
  common& cmn,
  int const& mlist)
{
  FEM_CMN_SVE(lulist);
  common_write write(cmn);
  int& n = cmn.n;
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_cref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_cref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<float> parf(cmn.parf, dimension(2000));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_cref<float> brat(cmn.brat, dimension(2000));
  arr_cref<int, 2> kfdp(cmn.kfdp, dimension(2000, 5));
  //
  str_arr_ref<1> chdl(sve.chdl, dimension(7));
  str_arr_ref<1> chmo(sve.chmo, dimension(12));
  if (is_called_first_time) {
    {
      static const char* values[] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
          "Sep", "Oct", "Nov", "Dec"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chmo;
    }
    {
      static const char* values[] = {
        "(())", " ", "()", "!!", "<>", "==", "(==)"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chdl;
    }
  }
  int lmx = fem::int0;
  int istr = fem::int0;
  int imax = fem::int0;
  int i = fem::int0;
  fem::str<16> chap = fem::char0;
  int len = fem::int0;
  int lem = fem::int0;
  int mdl = fem::int0;
  int ldl = fem::int0;
  fem::str<16> chac = fem::char0;
  int kc = fem::int0;
  int kcc = fem::int0;
  int j1 = fem::int0;
  int j2 = fem::int0;
  int j = fem::int0;
  int isep = fem::int0;
  arr_1d<6, float> ps(fem::fill0);
  int kf = fem::int0;
  fem::str<16> chan = fem::char0;
  int kfls = fem::int0;
  int kfla = fem::int0;
  int kflb = fem::int0;
  int kmul = fem::int0;
  int kflr = fem::int0;
  int kflc = fem::int0;
  int kflsp = fem::int0;
  int mstj24 = fem::int0;
  int kfmax = fem::int0;
  float pm = fem::float0;
  int idc = fem::int0;
  arr_1d<5, fem::str<16> > chad(fem::fill0);
  static const char* format_2700 = "(4x,i6,4x,a16,6x,i6,4x,a16)";
  //C
  //C...Purpose: to give program heading, or list an event, or particle
  //C...data, or current parameter values.
  //C
  //C...Initialization printout: version number and date of last change.
  //C      IF(MLIST.EQ.0.OR.MSTU(12).EQ.1) THEN
  //C        WRITE(MSTU(11),1000) MSTU(181),MSTU(182),MSTU(185),
  //C     &  CHMO(MSTU(184)),MSTU(183)
  //C        MSTU(12)=0
  //C        IF(MLIST.EQ.0) RETURN
  //C      ENDIF
  //C
  //C...List event data, including additional lines after N.
  if (mlist >= 1 && mlist <= 3) {
    if (mlist == 1) {
      write(mstu(11),
        "(/,/,/,28x,'Event listing (summary)',/,/,4x,'I  particle/jet KS',5x,"
        "'KF orig    p_x      p_y      p_z       E        m',/)");
    }
    if (mlist == 2) {
      write(mstu(11),
        "(/,/,/,28x,'Event listing (standard)',/,/,4x,'I  particle/jet',"
        "'  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',"
        "'       P(I,2)       P(I,3)       P(I,4)       P(I,5)',/)");
    }
    if (mlist == 3) {
      write(mstu(11),
        "(/,/,/,28x,'Event listing (with vertices)',/,/,4x,'I  particle/j',"
        "'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',"
        "'       P(I,2)       P(I,3)       P(I,4)       P(I,5)',/,73x,"
        "'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)',/)");
    }
    lmx = 12;
    if (mlist >= 2) {
      lmx = 16;
    }
    istr = 0;
    imax = n;
    if (mstu(2) > 0) {
      imax = mstu(2);
    }
    FEM_DO_SAFE(i, fem::max(1, mstu(1)), fem::max(imax, n + fem::max(0,
      mstu(3)))) {
      if ((i > imax && i <= n) || k(i, 1) < 0) {
        goto statement_120;
      }
      //C
      //C...Get particle name, pad it and check it is not too long.
      luname(cmn, k(i, 2), chap);
      len = 0;
      FEM_DO_SAFE(lem, 1, 16) {
        if (chap(lem, lem) != " ") {
          len = lem;
        }
      }
      mdl = (k(i, 1) + 19) / 10;
      ldl = 0;
      if (mdl == 2 || mdl >= 8) {
        chac = chap;
        if (len > lmx) {
          chac(lmx, lmx) = "?";
        }
      }
      else {
        ldl = 1;
        if (mdl == 1 || mdl == 7) {
          ldl = 2;
        }
        if (len == 0) {
          chac = chdl(mdl)(1, 2 * ldl) + str_cref(" ");
        }
        else {
          chac = chdl(mdl)(1, ldl) + chap(1, fem::min(len, lmx - 2 *
            ldl)) + chdl(mdl)(ldl + 1, 2 * ldl) + str_cref(" ");
          if (len + 2 * ldl > lmx) {
            chac(lmx, lmx) = "?";
          }
        }
      }
      //C
      //C...Add information on string connection.
      if (k(i, 1) == 1 || k(i, 1) == 2 || k(i, 1) == 11 || k(i, 1) == 12) {
        kc = lucomp(cmn, k(i, 2));
        kcc = 0;
        if (kc != 0) {
          kcc = kchg(kc, 2);
        }
        if (kcc != 0 && istr == 0) {
          istr = 1;
          if (len + 2 * ldl + 3 <= lmx) {
            chac(lmx - 1, lmx - 1) = "A";
          }
        }
        else if (kcc != 0 && (k(i, 1) == 2 || k(i, 1) == 12)) {
          if (len + 2 * ldl + 3 <= lmx) {
            chac(lmx - 1, lmx - 1) = "I";
          }
        }
        else if (kcc != 0) {
          istr = 0;
          if (len + 2 * ldl + 3 <= lmx) {
            chac(lmx - 1, lmx - 1) = "V";
          }
        }
      }
      //C
      //C...Write data for particle/jet.
      if (mlist == 1 && fem::abs(p(i, 4)) < 9999.f) {
        {
          write_loop wloop(cmn, mstu(11),
            "(1x,i4,2x,a12,1x,i2,1x,i6,1x,i4,5f9.3)");
          wloop, i, chac(1, 12);
          FEM_DO_SAFE(j1, 1, 3) {
            wloop, k(i, j1);
          }
          FEM_DO_SAFE(j2, 1, 5) {
            wloop, p(i, j2);
          }
        }
      }
      else if (mlist == 1 && fem::abs(p(i, 4)) < 99999.f) {
        {
          write_loop wloop(cmn, mstu(11),
            "(1x,i4,2x,a12,1x,i2,1x,i6,1x,i4,5f9.2)");
          wloop, i, chac(1, 12);
          FEM_DO_SAFE(j1, 1, 3) {
            wloop, k(i, j1);
          }
          FEM_DO_SAFE(j2, 1, 5) {
            wloop, p(i, j2);
          }
        }
      }
      else if (mlist == 1) {
        {
          write_loop wloop(cmn, mstu(11),
            "(1x,i4,2x,a12,1x,i2,1x,i6,1x,i4,5f9.1)");
          wloop, i, chac(1, 12);
          FEM_DO_SAFE(j1, 1, 3) {
            wloop, k(i, j1);
          }
          FEM_DO_SAFE(j2, 1, 5) {
            wloop, p(i, j2);
          }
        }
      }
      else if (mstu(5) == 10000 && (k(i, 1) == 3 || k(i, 1) == 13 || k(i,
        1) == 14)) {
        {
          write_loop wloop(cmn, mstu(11),
            "(1x,i4,2x,a16,1x,i3,1x,i8,2x,i4,2(3x,i1,2i4),5f13.5)");
          wloop, i, chac;
          FEM_DO_SAFE(j1, 1, 3) {
            wloop, k(i, j1);
          }
          wloop, k(i, 4) / 100000000, fem::mod(k(i, 4) / 10000,
            10000), fem::mod(k(i, 4), 10000), k(i, 5) / 100000000,
            fem::mod(k(i, 5) / 10000, 10000), fem::mod(k(i, 5),
            10000);
          FEM_DO_SAFE(j2, 1, 5) {
            wloop, p(i, j2);
          }
        }
      }
      else {
        {
          write_loop wloop(cmn, mstu(11),
            "(1x,i4,2x,a16,1x,i3,1x,i8,2x,i4,2(3x,i9),5f13.5)");
          wloop, i, chac;
          FEM_DO_SAFE(j1, 1, 5) {
            wloop, k(i, j1);
          }
          FEM_DO_SAFE(j2, 1, 5) {
            wloop, p(i, j2);
          }
        }
      }
      if (mlist == 3) {
        {
          write_loop wloop(cmn, mstu(11), "(66x,5(1x,f12.3))");
          FEM_DO_SAFE(j, 1, 5) {
            wloop, v(i, j);
          }
        }
      }
      //C
      //C...Insert extra separator lines specified by user.
      if (mstu(70) >= 1) {
        isep = 0;
        FEM_DO_SAFE(j, 1, fem::min(10, mstu(70))) {
          if (i == mstu(70 + j)) {
            isep = 1;
          }
        }
        if (isep == 1 && mlist == 1) {
          write(mstu(11), "(1x,78('='))");
        }
        if (isep == 1 && mlist >= 2) {
          write(mstu(11), "(1x,130('='))");
        }
      }
      statement_120:;
    }
    //C
    //C...Sum of charges and momenta.
    FEM_DO_SAFE(j, 1, 6) {
      ps(j) = plu(cmn, 0, j);
    }
    if (mlist == 1 && fem::abs(ps(4)) < 9999.f) {
      {
        write_loop wloop(cmn, mstu(11), "(19x,'sum:',f6.2,5x,5f9.3)");
        wloop, ps(6);
        FEM_DO_SAFE(j, 1, 5) {
          wloop, ps(j);
        }
      }
    }
    else if (mlist == 1 && fem::abs(ps(4)) < 99999.f) {
      {
        write_loop wloop(cmn, mstu(11), "(19x,'sum:',f6.2,5x,5f9.2)");
        wloop, ps(6);
        FEM_DO_SAFE(j, 1, 5) {
          wloop, ps(j);
        }
      }
    }
    else if (mlist == 1) {
      {
        write_loop wloop(cmn, mstu(11), "(19x,'sum:',f6.2,5x,5f9.1)");
        wloop, ps(6);
        FEM_DO_SAFE(j, 1, 5) {
          wloop, ps(j);
        }
      }
    }
    else {
      {
        write_loop wloop(cmn, mstu(11),
          "(19x,'sum charge:',f6.2,3x,'sum momentum and inv. mass:',5f13.5)");
        wloop, ps(6);
        FEM_DO_SAFE(j, 1, 5) {
          wloop, ps(j);
        }
      }
    }
    //C
    //C...Give simple list of KF codes defined in program.
  }
  else if (mlist == 11) {
    write(mstu(11), "(/,/,/,20x,'List of KF codes in program',/)");
    FEM_DO_SAFE(kf, 1, 40) {
      luname(cmn, kf, chap);
      luname(cmn, -kf, chan);
      if (chap != " " && chan == " ") {
        write(mstu(11), format_2700), kf, chap;
      }
      if (chan != " ") {
        write(mstu(11), format_2700), kf, chap, -kf, chan;
      }
    }
    FEM_DOSTEP(kfls, 1, 3, 2) {
      FEM_DO_SAFE(kfla, 1, 8) {
        FEM_DO_SAFE(kflb, 1, kfla - (3 - kfls) / 2) {
          kf = 1000 * kfla + 100 * kflb + kfls;
          luname(cmn, kf, chap);
          luname(cmn, -kf, chan);
          write(mstu(11), format_2700), kf, chap, -kf, chan;
        }
      }
    }
    FEM_DO_SAFE(kmul, 0, 5) {
      kfls = 3;
      if (kmul == 0 || kmul == 3) {
        kfls = 1;
      }
      if (kmul == 5) {
        kfls = 5;
      }
      kflr = 0;
      if (kmul == 2 || kmul == 3) {
        kflr = 1;
      }
      if (kmul == 4) {
        kflr = 2;
      }
      FEM_DO_SAFE(kflb, 1, 8) {
        FEM_DO_SAFE(kflc, 1, kflb - 1) {
          kf = 10000 * kflr + 100 * kflb + 10 * kflc + kfls;
          luname(cmn, kf, chap);
          luname(cmn, -kf, chan);
          write(mstu(11), format_2700), kf, chap, -kf, chan;
        }
        kf = 10000 * kflr + 110 * kflb + kfls;
        luname(cmn, kf, chap);
        write(mstu(11), format_2700), kf, chap;
      }
    }
    kf = 130;
    luname(cmn, kf, chap);
    write(mstu(11), format_2700), kf, chap;
    kf = 310;
    luname(cmn, kf, chap);
    write(mstu(11), format_2700), kf, chap;
    FEM_DO_SAFE(kflsp, 1, 3) {
      kfls = 2 + 2 * (kflsp / 3);
      FEM_DO_SAFE(kfla, 1, 8) {
        FEM_DO_SAFE(kflb, 1, kfla) {
          FEM_DO_SAFE(kflc, 1, kflb) {
            if (kflsp == 1 && (kfla == kflb || kflb == kflc)) {
              goto statement_180;
            }
            if (kflsp == 2 && kfla == kflc) {
              goto statement_180;
            }
            if (kflsp == 1) {
              kf = 1000 * kfla + 100 * kflc + 10 * kflb + kfls;
            }
            if (kflsp >= 2) {
              kf = 1000 * kfla + 100 * kflb + 10 * kflc + kfls;
            }
            luname(cmn, kf, chap);
            luname(cmn, -kf, chan);
            write(mstu(11), format_2700), kf, chap, -kf, chan;
            statement_180:;
          }
        }
      }
    }
    //C
    //C...List parton/particle data table. Check whether to be listed.
  }
  else if (mlist == 12) {
    write(mstu(11),
      "(/,/,/,30x,'Particle/parton data table',/,/,5x,'KF',5x,'KC',4x,"
      "'particle',8x,'antiparticle',6x,'chg  col  anti',8x,'mass',7x,'width',"
      "7x,'w-cut',5x,'lifetime',1x,'decay',/,11x,'IDC',1x,'on/off',1x,'ME',3x,"
      "'Br.rat.',4x,'decay products')");
    mstj24 = mstj(24);
    mstj(24) = 0;
    kfmax = 20883;
    if (mstu(2) != 0) {
      kfmax = mstu(2);
    }
    FEM_DO_SAFE(kf, fem::max(1, mstu(1)), kfmax) {
      kc = lucomp(cmn, kf);
      if (kc == 0) {
        goto statement_220;
      }
      if (mstu(14) == 0 && kf > 100 && kc <= 100) {
        goto statement_220;
      }
      if (mstu(14) > 0 && kf > 100 && fem::max(fem::mod(kf / 1000,
          10), fem::mod(kf / 100, 10)) > mstu(14)) {
        goto statement_220;
      }
      //C
      //C...Find particle name and mass. Print information.
      luname(cmn, kf, chap);
      if (kf <= 100 && chap == " " && mdcy(kc, 2) == 0) {
        goto statement_220;
      }
      luname(cmn, -kf, chan);
      pm = ulmass(cmn, kf);
      write(mstu(11),
        "(/,1x,i6,3x,i4,4x,a16,a16,3i5,1x,f12.5,2(1x,f11.5),2x,f12.5,3x,i2)"),
        kf, kc, chap, chan, kchg(kc, 1), kchg(kc, 2), kchg(kc, 3),
        pm, pmas(kc, 2), pmas(kc, 3), pmas(kc, 4), mdcy(kc, 1);
      //C
      //C...Particle decay: channel number, branching ration, matrix element,
      //C...decay products.
      if (kf > 100 && kc <= 100) {
        goto statement_220;
      }
      FEM_DO_SAFE(idc, mdcy(kc, 2), mdcy(kc, 2) + mdcy(kc, 3) - 1) {
        FEM_DO_SAFE(j, 1, 5) {
          luname(cmn, kfdp(idc, j), chad(j));
        }
        {
          write_loop wloop(cmn, mstu(11),
            "(10x,i4,2x,i3,2x,i3,2x,f8.5,4x,5a16)");
          wloop, idc, mdme(idc, 1), mdme(idc, 2), brat(idc);
          FEM_DO_SAFE(j, 1, 5) {
            wloop, chad(j);
          }
        }
      }
      statement_220:;
    }
    mstj(24) = mstj24;
    //C
    //C...List parameter value table.
  }
  else if (mlist == 13) {
    write(mstu(11),
      "(/,/,/,20x,'Parameter value table',/,/,4x,'I',3x,'MSTU(I)',8x,"
      "'PARU(I)',3x,'MSTJ(I)',8x,'PARJ(I)',8x,'PARF(I)')");
    FEM_DO_SAFE(i, 1, 200) {
      write(mstu(11),
        "(1x,i4,1x,i9,1x,f14.5,1x,i9,1x,f14.5,1x,f14.5)"), i, mstu(i),
        paru(i), mstj(i), parj(i), parf(i);
    }
  }
  //C
  //C...Format statements for output on unit MSTU(11) (by default 6).
  //Clin 1000 FORMAT(///20X,'The Lund Monte Carlo - JETSET version ',I1,'.',I1/
  //Clin     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)
  //C
}

//C
//C*********************************************************************
//C
void
luerrm(
  common& cmn,
  int const& merr,
  str_cref chmess)
{
  common_write write(cmn);
  // COMMON ludat1
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  //
  static const char* format_1100 =
    "(/,5x,'Error type',i2,' has occured after',i6,' LUEXEC calls:',/,5x,a)";
  //C
  //C...Purpose: to inform user of errors in program execution.
  //C
  write(6, star), "merr,chmess=", merr, chmess;
  //C
  //C...Write first few warnings, then be silent.
  if (merr <= 10) {
    mstu(27)++;
    mstu(28) = merr;
    if (mstu(25) == 1 && mstu(27) <= mstu(26)) {
      write(mstu(11),
        "(/,5x,'Advisory warning type',i2,' given after',i6,' LUEXEC calls:',"
        "/,5x,a)"),
        merr, mstu(31), chmess;
    }
    //C
    //C...Write first few errors, then be silent or stop program.
  }
  else if (merr <= 20) {
    mstu(23)++;
    mstu(24) = merr - 10;
    if (mstu(21) >= 1 && mstu(23) <= mstu(22)) {
      write(mstu(11), format_1100), merr - 10, mstu(31), chmess;
    }
    if (mstu(21) >= 2 && mstu(23) > mstu(22)) {
      write(mstu(11), format_1100), merr - 10, mstu(31), chmess;
      write(mstu(11),
        "(5x,'Execution will be stopped after listing of last ','event!')");
      if (merr != 17) {
        lulist(cmn, 2);
      }
      FEM_STOP(0);
    }
    //C
    //C...Stop program in case of irreparable error.
  }
  else {
    write(mstu(11),
      "(/,5x,'Fatal error type',i2,' has occured after',i6,' LUEXEC calls:',/,"
      "5x,a,/,5x,'Execution will now be stopped!')"),
      merr - 20, mstu(31), chmess;
    FEM_STOP(0);
  }
  //C
  //C...Formats for output.
  //C
}

//C
//C*********************************************************************
//C
void
lukfdi(
  common& cmn,
  int const& kfl1,
  int const& kfl2,
  int& kfl3,
  int& kf)
{
  arr_cref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<float> parf(cmn.parf, dimension(2000));
  //
  int kf1a = fem::int0;
  int kf2a = fem::int0;
  int ktab1 = fem::int0;
  int kfl1a = fem::int0;
  int kfl1b = fem::int0;
  int kfl1s = fem::int0;
  int ktab2 = fem::int0;
  int kfl2a = fem::int0;
  int kfl2b = fem::int0;
  int kfl2s = fem::int0;
  float par2 = fem::float0;
  float par3 = fem::float0;
  float par4 = fem::float0;
  float par3m = fem::float0;
  float par4m = fem::float0;
  float pardm = fem::float0;
  float pars0 = fem::float0;
  float pars1 = fem::float0;
  float pars2 = fem::float0;
  float parsm = fem::float0;
  int mbary = fem::int0;
  int kfda = fem::int0;
  int kflda = fem::int0;
  int kfldb = fem::int0;
  int kflds = fem::int0;
  float wtdq = fem::float0;
  int kfs = fem::int0;
  int kfla = fem::int0;
  int kflb = fem::int0;
  int kfl1d = fem::int0;
  int kfl1e = fem::int0;
  int kfl3a = fem::int0;
  int kmul = fem::int0;
  float rmul = fem::float0;
  int kfls = fem::int0;
  float rmix = fem::float0;
  int imix = fem::int0;
  int kflc = fem::int0;
  int kbary = fem::int0;
  float wt = fem::float0;
  int kfld = fem::int0;
  int kflf = fem::int0;
  int kfle = fem::int0;
  int kfll = fem::int0;
  int kt3l = fem::int0;
  int kt3u = fem::int0;
  float rfl = fem::float0;
  int kts = fem::int0;
  int kt3 = fem::int0;
  int ktabs = fem::int0;
  int ktab3 = fem::int0;
  int kfl3b = fem::int0;
  int kc = fem::int0;
  //C
  //C...Purpose: to generate a new flavour pair and combine off a hadron.
  //C
  //C...Default flavour values. Input consistency checks.
  kf1a = fem::iabs(kfl1);
  kf2a = fem::iabs(kfl2);
  kfl3 = 0;
  kf = 0;
  if (kf1a == 0) {
    return;
  }
  if (kf2a != 0) {
    if (kf1a <= 10 && kf2a <= 10 && kfl1 * kfl2 > 0) {
      return;
    }
    if (kf1a > 10 && kf2a > 10) {
      return;
    }
    if ((kf1a > 10 || kf2a > 10) && kfl1 * kfl2 < 0) {
      return;
    }
  }
  //C
  //C...Check if tabulated flavour probabilities are to be used.
  if (mstj(15) == 1) {
    ktab1 = -1;
    if (kf1a >= 1 && kf1a <= 6) {
      ktab1 = kf1a;
    }
    kfl1a = fem::mod(kf1a / 1000, 10);
    kfl1b = fem::mod(kf1a / 100, 10);
    kfl1s = fem::mod(kf1a, 10);
    if (kfl1a >= 1 && kfl1a <= 4 && kfl1b >= 1 && kfl1b <= 4) {
      ktab1 = 6 + kfl1a * (kfl1a - 2) + 2 * kfl1b + (kfl1s - 1) / 2;
    }
    if (kfl1a >= 1 && kfl1a <= 4 && kfl1a == kfl1b) {
      ktab1 = ktab1 - 1;
    }
    if (kf1a >= 1 && kf1a <= 6) {
      kfl1a = kf1a;
    }
    ktab2 = 0;
    if (kf2a != 0) {
      ktab2 = -1;
      if (kf2a >= 1 && kf2a <= 6) {
        ktab2 = kf2a;
      }
      kfl2a = fem::mod(kf2a / 1000, 10);
      kfl2b = fem::mod(kf2a / 100, 10);
      kfl2s = fem::mod(kf2a, 10);
      if (kfl2a >= 1 && kfl2a <= 4 && kfl2b >= 1 && kfl2b <= 4) {
        ktab2 = 6 + kfl2a * (kfl2a - 2) + 2 * kfl2b + (kfl2s - 1) / 2;
      }
      if (kfl2a >= 1 && kfl2a <= 4 && kfl2a == kfl2b) {
        ktab2 = ktab2 - 1;
      }
    }
    if (ktab1 >= 0 && ktab2 >= 0) {
      goto statement_140;
    }
  }
  //C
  //C...Parameters and breaking diquark parameter combinations.
  statement_100:
  par2 = parj(2);
  par3 = parj(3);
  par4 = 3.f * parj(4);
  if (mstj(12) >= 2) {
    par3m = fem::sqrt(parj(3));
    par4m = 1.f / (3.f * fem::sqrt(parj(4)));
    pardm = parj(7) / (parj(7) + par3m * parj(6));
    pars0 = parj(5) * (2.f + (1.f + par2 * par3m * parj(7)) * (1.f + par4m));
    pars1 = parj(7) * pars0 / (2.f * par3m) + parj(5) * (parj(6) * (
      1.f + par4m) + par2 * par3m * parj(6) * parj(7));
    pars2 = parj(5) * 2.f * parj(6) * parj(7) * (par2 * parj(7) + (
      1.f + par4m) / par3m);
    parsm = fem::max(pars0, pars1, pars2);
    par4 = par4 * (1.f + parsm) / (1.f + parsm / (3.f * par4m));
  }
  //C
  //C...Choice of whether to generate meson or baryon.
  mbary = 0;
  kfda = 0;
  if (kf1a <= 10) {
    if (kf2a == 0 && mstj(12) >= 1 && (1.f + parj(1)) * rlu(cmn, 0) > 1.f) {
      mbary = 1;
    }
    if (kf2a > 10) {
      mbary = 2;
    }
    if (kf2a > 10 && kf2a <= 10000) {
      kfda = kf2a;
    }
  }
  else {
    mbary = 2;
    if (kf1a <= 10000) {
      kfda = kf1a;
    }
  }
  //C
  //C...Possibility of process diquark -> meson + new diquark.
  if (kfda != 0 && mstj(12) >= 2) {
    kflda = fem::mod(kfda / 1000, 10);
    kfldb = fem::mod(kfda / 100, 10);
    kflds = fem::mod(kfda, 10);
    wtdq = pars0;
    if (fem::max(kflda, kfldb) == 3) {
      wtdq = pars1;
    }
    if (fem::min(kflda, kfldb) == 3) {
      wtdq = pars2;
    }
    if (kflds == 1) {
      wtdq = wtdq / (3.f * par4m);
    }
    if ((1.f + wtdq) * rlu(cmn, 0) > 1.f) {
      mbary = -1;
    }
    if (mbary ==  - 1 && kf2a != 0) {
      return;
    }
  }
  //C
  //C...Flavour for meson, possibly with new flavour.
  if (mbary <= 0) {
    kfs = fem::isign(1, kfl1);
    if (mbary == 0) {
      if (kf2a == 0) {
        kfl3 = fem::isign(1 + fem::fint((2.f + par2) * rlu(cmn, 0)), -kfl1);
      }
      kfla = fem::max(kf1a, kf2a + fem::iabs(kfl3));
      kflb = fem::min(kf1a, kf2a + fem::iabs(kfl3));
      if (kfla != kf1a) {
        kfs = -kfs;
      }
      //C
      //C...Splitting of diquark into meson plus new diquark.
    }
    else {
      kfl1a = fem::mod(kf1a / 1000, 10);
      kfl1b = fem::mod(kf1a / 100, 10);
      statement_110:
      kfl1d = kfl1a + fem::fint(rlu(cmn, 0) + 0.5f) * (kfl1b - kfl1a);
      kfl1e = kfl1a + kfl1b - kfl1d;
      if ((kfl1d == 3 && rlu(cmn, 0) > pardm) || (kfl1e == 3 && rlu(cmn,
          0) < pardm)) {
        kfl1d = kfl1a + kfl1b - kfl1d;
        kfl1e = kfl1a + kfl1b - kfl1e;
      }
      kfl3a = 1 + fem::fint((2.f + par2 * par3m * parj(7)) * rlu(cmn, 0));
      if ((kfl1e != kfl3a && rlu(cmn, 0) > (1.f + par4m) / fem::max(2.f,
          1.f + par4m)) || (kfl1e == kfl3a && rlu(cmn, 0) > 2.f / fem::max(2.f,
          1.f + par4m))) {
        goto statement_110;
      }
      kflds = 3;
      if (kfl1e != kfl3a) {
        kflds = 2 * fem::fint(rlu(cmn, 0) + 1.f / (1.f + par4m)) + 1;
      }
      kfl3 = fem::isign(10000 + 1000 * fem::max(kfl1e, kfl3a) + 100 *
        fem::min(kfl1e, kfl3a) + kflds, -kfl1);
      kfla = fem::max(kfl1d, kfl3a);
      kflb = fem::min(kfl1d, kfl3a);
      if (kfla != kfl1d) {
        kfs = -kfs;
      }
    }
    //C
    //C...Form meson, with spin and flavour mixing for diagonal states.
    if (kfla <= 2) {
      kmul = fem::fint(parj(11) + rlu(cmn, 0));
    }
    if (kfla == 3) {
      kmul = fem::fint(parj(12) + rlu(cmn, 0));
    }
    if (kfla >= 4) {
      kmul = fem::fint(parj(13) + rlu(cmn, 0));
    }
    if (kmul == 0 && parj(14) > 0.f) {
      if (rlu(cmn, 0) < parj(14)) {
        kmul = 2;
      }
    }
    else if (kmul == 1 && parj(15) + parj(16) + parj(17) > 0.f) {
      rmul = rlu(cmn, 0);
      if (rmul < parj(15)) {
        kmul = 3;
      }
      if (kmul == 1 && rmul < parj(15) + parj(16)) {
        kmul = 4;
      }
      if (kmul == 1 && rmul < parj(15) + parj(16) + parj(17)) {
        kmul = 5;
      }
    }
    kfls = 3;
    if (kmul == 0 || kmul == 3) {
      kfls = 1;
    }
    if (kmul == 5) {
      kfls = 5;
    }
    if (kfla != kflb) {
      kf = (100 * kfla + 10 * kflb + kfls) * kfs * fem::pow((-1), kfla);
    }
    else {
      rmix = rlu(cmn, 0);
      imix = 2 * kfla + 10 * kmul;
      if (kfla <= 3) {
        kf = 110 * (1 + fem::fint(rmix + parf(imix - 1)) + fem::fint(
          rmix + parf(imix))) + kfls;
      }
      if (kfla >= 4) {
        kf = 110 * kfla + kfls;
      }
    }
    if (kmul == 2 || kmul == 3) {
      kf += fem::isign(10000, kf);
    }
    if (kmul == 4) {
      kf += fem::isign(20000, kf);
    }
    //C
    //C...Generate diquark flavour.
  }
  else {
    statement_120:
    if (kf1a <= 10 && kf2a == 0) {
      kfla = kf1a;
      statement_130:
      kflb = 1 + fem::fint((2.f + par2 * par3) * rlu(cmn, 0));
      kflc = 1 + fem::fint((2.f + par2 * par3) * rlu(cmn, 0));
      kflds = 1;
      if (kflb >= kflc) {
        kflds = 3;
      }
      if (kflds == 1 && par4 * rlu(cmn, 0) > 1.f) {
        goto statement_130;
      }
      if (kflds == 3 && par4 < rlu(cmn, 0)) {
        goto statement_130;
      }
      kfl3 = fem::isign(1000 * fem::max(kflb, kflc) + 100 * fem::min(kflb,
        kflc) + kflds, kfl1);
      //C
      //C...Take diquark flavour from input.
    }
    else if (kf1a <= 10) {
      kfla = kf1a;
      kflb = fem::mod(kf2a / 1000, 10);
      kflc = fem::mod(kf2a / 100, 10);
      kflds = fem::mod(kf2a, 10);
      //C
      //C...Generate (or take from input) quark to go with diquark.
    }
    else {
      if (kf2a == 0) {
        kfl3 = fem::isign(1 + fem::fint((2.f + par2) * rlu(cmn, 0)), kfl1);
      }
      kfla = kf2a + fem::iabs(kfl3);
      kflb = fem::mod(kf1a / 1000, 10);
      kflc = fem::mod(kf1a / 100, 10);
      kflds = fem::mod(kf1a, 10);
    }
    //C
    //C...SU(6) factors for formation of baryon. Try again if fails.
    kbary = kflds;
    if (kflds == 3 && kflb != kflc) {
      kbary = 5;
    }
    if (kfla != kflb && kfla != kflc) {
      kbary++;
    }
    wt = parf(60 + kbary) + parj(18) * parf(70 + kbary);
    if (mbary == 1 && mstj(12) >= 2) {
      wtdq = pars0;
      if (fem::max(kflb, kflc) == 3) {
        wtdq = pars1;
      }
      if (fem::min(kflb, kflc) == 3) {
        wtdq = pars2;
      }
      if (kflds == 1) {
        wtdq = wtdq / (3.f * par4m);
      }
      if (kflds == 1) {
        wt = wt * (1.f + wtdq) / (1.f + parsm / (3.f * par4m));
      }
      if (kflds == 3) {
        wt = wt * (1.f + wtdq) / (1.f + parsm);
      }
    }
    if (kf2a == 0 && wt < rlu(cmn, 0)) {
      goto statement_120;
    }
    //C
    //C...Form baryon. Distinguish Lambda- and Sigmalike baryons.
    kfld = fem::max(kfla, kflb, kflc);
    kflf = fem::min(kfla, kflb, kflc);
    kfle = kfla + kflb + kflc - kfld - kflf;
    kfls = 2;
    if ((parf(60 + kbary) + parj(18) * parf(70 + kbary)) * rlu(cmn,
        0) > parf(60 + kbary)) {
      kfls = 4;
    }
    kfll = 0;
    if (kfls == 2 && kfld > kfle && kfle > kflf) {
      if (kflds == 1 && kfla == kfld) {
        kfll = 1;
      }
      if (kflds == 1 && kfla != kfld) {
        kfll = fem::fint(0.25f + rlu(cmn, 0));
      }
      if (kflds == 3 && kfla != kfld) {
        kfll = fem::fint(0.75f + rlu(cmn, 0));
      }
    }
    if (kfll == 0) {
      kf = fem::isign(1000 * kfld + 100 * kfle + 10 * kflf + kfls, kfl1);
    }
    if (kfll == 1) {
      kf = fem::isign(1000 * kfld + 100 * kflf + 10 * kfle + kfls, kfl1);
    }
  }
  return;
  //C
  //C...Use tabulated probabilities to select new flavour and hadron.
  statement_140:
  if (ktab2 == 0 && mstj(12) <= 0) {
    kt3l = 1;
    kt3u = 6;
  }
  else if (ktab2 == 0 && ktab1 >= 7 && mstj(12) <= 1) {
    kt3l = 1;
    kt3u = 6;
  }
  else if (ktab2 == 0) {
    kt3l = 1;
    kt3u = 22;
  }
  else {
    kt3l = ktab2;
    kt3u = ktab2;
  }
  rfl = 0.f;
  FEM_DO_SAFE(kts, 0, 2) {
    FEM_DO_SAFE(kt3, kt3l, kt3u) {
      rfl += parf(120 + 80 * ktab1 + 25 * kts + kt3);
    }
  }
  rfl = rlu(cmn, 0) * rfl;
  FEM_DO_SAFE(kts, 0, 2) {
    ktabs = kts;
    FEM_DO_SAFE(kt3, kt3l, kt3u) {
      ktab3 = kt3;
      rfl = rfl - parf(120 + 80 * ktab1 + 25 * kts + kt3);
      if (rfl <= 0.f) {
        goto statement_170;
      }
    }
  }
  statement_170:
  //C
  //C...Reconstruct flavour of produced quark/diquark.
  if (ktab3 <= 6) {
    kfl3a = ktab3;
    kfl3b = 0;
    kfl3 = fem::isign(kfl3a, kfl1 * (2 * ktab1 - 13));
  }
  else {
    kfl3a = 1;
    if (ktab3 >= 8) {
      kfl3a = 2;
    }
    if (ktab3 >= 11) {
      kfl3a = 3;
    }
    if (ktab3 >= 16) {
      kfl3a = 4;
    }
    kfl3b = (ktab3 - 6 - kfl3a * (kfl3a - 2)) / 2;
    kfl3 = 1000 * kfl3a + 100 * kfl3b + 1;
    if (kfl3a == kfl3b || ktab3 != 6 + kfl3a * (kfl3a - 2) + 2 * kfl3b) {
      kfl3 += 2;
    }
    kfl3 = fem::isign(kfl3, kfl1 * (13 - 2 * ktab1));
  }
  //C
  //C...Reconstruct meson code.
  if (kfl3a == kfl1a && kfl3b == kfl1b && (kfl3a <= 3 || kfl3b != 0)) {
    rfl = rlu(cmn, 0) * (parf(143 + 80 * ktab1 + 25 * ktabs) + parf(
      144 + 80 * ktab1 + 25 * ktabs) + parf(145 + 80 * ktab1 + 25 *
      ktabs));
    kf = 110 + 2 * ktabs + 1;
    if (rfl > parf(143 + 80 * ktab1 + 25 * ktabs)) {
      kf = 220 + 2 * ktabs + 1;
    }
    if (rfl > parf(143 + 80 * ktab1 + 25 * ktabs) + parf(144 + 80 *
        ktab1 + 25 * ktabs)) {
      kf = 330 + 2 * ktabs + 1;
    }
  }
  else if (ktab1 <= 6 && ktab3 <= 6) {
    kfla = fem::max(ktab1, ktab3);
    kflb = fem::min(ktab1, ktab3);
    kfs = fem::isign(1, kfl1);
    if (kfla != kf1a) {
      kfs = -kfs;
    }
    kf = (100 * kfla + 10 * kflb + 2 * ktabs + 1) * kfs * fem::pow((-1), kfla);
  }
  else if (ktab1 >= 7 && ktab3 >= 7) {
    kfs = fem::isign(1, kfl1);
    if (kfl1a == kfl3a) {
      kfla = fem::max(kfl1b, kfl3b);
      kflb = fem::min(kfl1b, kfl3b);
      if (kfla != kfl1b) {
        kfs = -kfs;
      }
    }
    else if (kfl1a == kfl3b) {
      kfla = kfl3a;
      kflb = kfl1b;
      kfs = -kfs;
    }
    else if (kfl1b == kfl3a) {
      kfla = kfl1a;
      kflb = kfl3b;
    }
    else if (kfl1b == kfl3b) {
      kfla = fem::max(kfl1a, kfl3a);
      kflb = fem::min(kfl1a, kfl3a);
      if (kfla != kfl1a) {
        kfs = -kfs;
      }
    }
    else {
      luerrm(cmn, 2, "(LUKFDI:) no matching flavours for qq -> qq");
      goto statement_100;
    }
    kf = (100 * kfla + 10 * kflb + 2 * ktabs + 1) * kfs * fem::pow((-1), kfla);
    //C
    //C...Reconstruct baryon code.
  }
  else {
    if (ktab1 >= 7) {
      kfla = kfl3a;
      kflb = kfl1a;
      kflc = kfl1b;
    }
    else {
      kfla = kfl1a;
      kflb = kfl3a;
      kflc = kfl3b;
    }
    kfld = fem::max(kfla, kflb, kflc);
    kflf = fem::min(kfla, kflb, kflc);
    kfle = kfla + kflb + kflc - kfld - kflf;
    if (ktabs == 0) {
      kf = fem::isign(1000 * kfld + 100 * kflf + 10 * kfle + 2, kfl1);
    }
    if (ktabs >= 1) {
      kf = fem::isign(1000 * kfld + 100 * kfle + 10 * kflf + 2 * ktabs, kfl1);
    }
  }
  //C
  //C...Check that constructed flavour code is an allowed one.
  if (kfl2 != 0) {
    kfl3 = 0;
  }
  kc = lucomp(cmn, kf);
  if (kc == 0) {
    luerrm(cmn, 2, "(LUKFDI:) user-defined flavour probabilities " +
      str_cref("failed"));
    goto statement_100;
  }
  //C
}

//C
//C*********************************************************************
//C
void
luprep(
  common& cmn,
  int const& ip)
{
  common_write write(cmn);
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  int i1 = fem::int0;
  int mqgst = fem::int0;
  int i = fem::int0;
  int kc = fem::int0;
  int kq = fem::int0;
  int kcs = fem::int0;
  int ia = fem::int0;
  int nstp = fem::int0;
  int j = fem::int0;
  int ib = fem::int0;
  int mrev = fem::int0;
  int ns = fem::int0;
  int nsin = fem::int0;
  float pdm = fem::float0;
  int ic = fem::int0;
  arr_1d<5, double> dps(fem::fill0);
  float pd = fem::float0;
  arr_1d<5, double> dpc(fem::fill0);
  int ic1 = fem::int0;
  int ic2 = fem::int0;
  int nsav = fem::int0;
  float pecm = fem::float0;
  int kc1 = fem::int0;
  int kc2 = fem::int0;
  int kq1 = fem::int0;
  int kq2 = fem::int0;
  int kfln = fem::int0;
  int kfldmp = fem::int0;
  int kfdmp = fem::int0;
  int kflm = fem::int0;
  float pa = fem::float0;
  arr_1d<3, float> ue(fem::fill0);
  float phi = fem::float0;
  int np = fem::int0;
  float ha = fem::float0;
  float hd1 = fem::float0;
  float hd2 = fem::float0;
  float hr = fem::float0;
  float hc = fem::float0;
  float hk1 = fem::float0;
  float hk2 = fem::float0;
  int ir = fem::int0;
  int mcomb = fem::int0;
  int kci = fem::int0;
  float hcr = fem::float0;
  float hb = fem::float0;
  float hd = fem::float0;
  int kfn = fem::int0;
  int kqs = fem::int0;
  //C
  //C...Purpose: to rearrange partons along strings, to allow small systems
  //C...to collapse into one or two particles and to check flavours.
  //C
  //C...Rearrange parton shower product listing along strings: begin loop.
  i1 = n;
  FEM_DO_SAFE(mqgst, 1, 2) {
    FEM_DO_SAFE(i, fem::max(1, ip), n) {
      if (k(i, 1) != 3) {
        goto statement_120;
      }
      kc = lucomp(cmn, k(i, 2));
      if (kc == 0) {
        goto statement_120;
      }
      kq = kchg(kc, 2);
      if (kq == 0 || (mqgst == 1 && kq == 2)) {
        goto statement_120;
      }
      //C
      //C...Pick up loose string end.
      kcs = 4;
      if (kq * fem::isign(1, k(i, 2)) < 0) {
        kcs = 5;
      }
      ia = i;
      nstp = 0;
      statement_100:
      nstp++;
      if (nstp > 4 * n) {
        luerrm(cmn, 14, "(LUPREP:) caught in infinite loop");
        return;
      }
      //C
      //C...Copy undecayed parton.
      if (k(ia, 1) == 3) {
        if (i1 >= mstu(4) - mstu(32) - 5) {
          luerrm(cmn, 11, "(LUPREP:) no more memory left in LUJETS");
          return;
        }
        i1++;
        k(i1, 1) = 2;
        if (nstp >= 2 && fem::iabs(k(ia, 2)) != 21) {
          k(i1, 1) = 1;
        }
        k(i1, 2) = k(ia, 2);
        k(i1, 3) = ia;
        k(i1, 4) = 0;
        k(i1, 5) = 0;
        FEM_DO_SAFE(j, 1, 5) {
          p(i1, j) = p(ia, j);
          v(i1, j) = v(ia, j);
        }
        k(ia, 1) += 10;
        if (k(i1, 1) == 1) {
          goto statement_120;
        }
      }
      //C
      //C...Go to next parton in colour space.
      ib = ia;
      if (fem::mod(k(ib, kcs) / fem::pow2(mstu(5)), 2) == 0 && fem::mod(k(ib,
          kcs), mstu(5)) != 0) {
        ia = fem::mod(k(ib, kcs), mstu(5));
        k(ib, kcs) += fem::pow2(mstu(5));
        mrev = 0;
      }
      else {
        if (k(ib, kcs) >= 2 * fem::pow2(mstu(5)) || fem::mod(k(ib,
            kcs) / mstu(5), mstu(5)) == 0) {
          kcs = 9 - kcs;
        }
        ia = fem::mod(k(ib, kcs) / mstu(5), mstu(5));
        k(ib, kcs) += 2 * fem::pow2(mstu(5));
        mrev = 1;
      }
      if (ia <= 0 || ia > n) {
        luerrm(cmn, 12, "(LUPREP:) colour rearrangement failed");
        return;
      }
      if (fem::mod(k(ia, 4) / mstu(5), mstu(5)) == ib || fem::mod(k(ia,
          5) / mstu(5), mstu(5)) == ib) {
        if (mrev == 1) {
          kcs = 9 - kcs;
        }
        if (fem::mod(k(ia, kcs) / mstu(5), mstu(5)) != ib) {
          kcs = 9 - kcs;
        }
        k(ia, kcs) += 2 * fem::pow2(mstu(5));
      }
      else {
        if (mrev == 0) {
          kcs = 9 - kcs;
        }
        if (fem::mod(k(ia, kcs), mstu(5)) != ib) {
          kcs = 9 - kcs;
        }
        k(ia, kcs) += fem::pow2(mstu(5));
      }
      if (ia != i) {
        goto statement_100;
      }
      k(i1, 1) = 1;
      statement_120:;
    }
  }
  n = i1;
  //C
  //C...Find lowest-mass colour singlet jet system, OK if above thresh.
  if (mstj(14) <= 0) {
    goto statement_320;
  }
  ns = n;
  statement_140:
  nsin = n - ns;
  pdm = 1.f + parj(32);
  ic = 0;
  FEM_DO_SAFE(i, fem::max(1, ip), ns) {
    if (k(i, 1) != 1 && k(i, 1) != 2) {
    }
    else if (k(i, 1) == 2 && ic == 0) {
      nsin++;
      ic = i;
      FEM_DO_SAFE(j, 1, 4) {
        dps(j) = fem::dble(p(i, j));
      }
      mstj(93) = 1;
      dps(5) = fem::dble(ulmass(cmn, k(i, 2)));
    }
    else if (k(i, 1) == 2) {
      FEM_DO_SAFE(j, 1, 4) {
        dps(j) += fem::dble(p(i, j));
      }
    }
    else if (ic != 0 && kchg(lucomp(cmn, k(i, 2)), 2) != 0) {
      FEM_DO_SAFE(j, 1, 4) {
        dps(j) += fem::dble(p(i, j));
      }
      mstj(93) = 1;
      dps(5) += fem::dble(ulmass(cmn, k(i, 2)));
      pd = fem::sngl(fem::sqrt(fem::max(0e0, fem::pow2(dps(4)) -
        fem::pow2(dps(1)) - fem::pow2(dps(2)) - fem::pow2(dps(3)))) -
        dps(5));
      if (pd < pdm) {
        pdm = pd;
        FEM_DO_SAFE(j, 1, 5) {
          dpc(j) = dps(j);
        }
        ic1 = ic;
        ic2 = i;
      }
      ic = 0;
    }
    else {
      nsin++;
    }
  }
  if (pdm >= parj(32)) {
    goto statement_320;
  }
  //C
  //C...Fill small-mass system as cluster.
  nsav = n;
  pecm = fem::sngl(fem::sqrt(fem::max(0e0, fem::pow2(dpc(4)) -
    fem::pow2(dpc(1)) - fem::pow2(dpc(2)) - fem::pow2(dpc(3)))));
  k(n + 1, 1) = 11;
  k(n + 1, 2) = 91;
  k(n + 1, 3) = ic1;
  k(n + 1, 4) = n + 2;
  k(n + 1, 5) = n + 3;
  p(n + 1, 1) = fem::sngl(dpc(1));
  p(n + 1, 2) = fem::sngl(dpc(2));
  p(n + 1, 3) = fem::sngl(dpc(3));
  p(n + 1, 4) = fem::sngl(dpc(4));
  p(n + 1, 5) = pecm;
  //C
  //C...Form two particles from flavours of lowest-mass system, if feasible.
  k(n + 2, 1) = 1;
  k(n + 3, 1) = 1;
  if (mstu(16) != 2) {
    k(n + 2, 3) = n + 1;
    k(n + 3, 3) = n + 1;
  }
  else {
    k(n + 2, 3) = ic1;
    k(n + 3, 3) = ic2;
  }
  k(n + 2, 4) = 0;
  k(n + 3, 4) = 0;
  k(n + 2, 5) = 0;
  k(n + 3, 5) = 0;
  if (fem::iabs(k(ic1, 2)) != 21) {
    kc1 = lucomp(cmn, k(ic1, 2));
    kc2 = lucomp(cmn, k(ic2, 2));
    if (kc1 == 0 || kc2 == 0) {
      goto statement_320;
    }
    kq1 = kchg(kc1, 2) * fem::isign(1, k(ic1, 2));
    kq2 = kchg(kc2, 2) * fem::isign(1, k(ic2, 2));
    if (kq1 + kq2 != 0) {
      goto statement_320;
    }
    statement_200:
    lukfdi(cmn, k(ic1, 2), 0, kfln, k(n + 2, 2));
    lukfdi(cmn, k(ic2, 2), -kfln, kfldmp, k(n + 3, 2));
    if (k(n + 2, 2) == 0 || k(n + 3, 2) == 0) {
      goto statement_200;
    }
  }
  else {
    if (fem::iabs(k(ic2, 2)) != 21) {
      goto statement_320;
    }
    statement_210:
    lukfdi(cmn, 1 + fem::fint((2.f + parj(2)) * rlu(cmn, 0)), 0, kfln, kfdmp);
    lukfdi(cmn, kfln, 0, kflm, k(n + 2, 2));
    lukfdi(cmn, -kfln, -kflm, kfldmp, k(n + 3, 2));
    if (k(n + 2, 2) == 0 || k(n + 3, 2) == 0) {
      goto statement_210;
    }
  }
  p(n + 2, 5) = ulmass(cmn, k(n + 2, 2));
  p(n + 3, 5) = ulmass(cmn, k(n + 3, 2));
  if (p(n + 2, 5) + p(n + 3, 5) + parj(64) >= pecm && nsin == 1) {
    goto statement_320;
  }
  if (p(n + 2, 5) + p(n + 3, 5) + parj(64) >= pecm) {
    goto statement_260;
  }
  //C
  //C...Perform two-particle decay of jet system, if possible.
  //Clin-5/2012:
  //C      IF(PECM.GE.0.02d0*DPC(4)) THEN
  if (fem::dble(pecm) >= 0.02e0 * dpc(4)) {
    pa = fem::sqrt((fem::pow2(pecm) - fem::pow2((p(n + 2, 5) + p(n + 3,
      5)))) * (fem::pow2(pecm) - fem::pow2((p(n + 2, 5) - p(n + 3,
      5))))) / (2.f * pecm);
    ue(3) = 2.f * rlu(cmn, 0) - 1.f;
    phi = paru(2) * rlu(cmn, 0);
    ue(1) = fem::sqrt(1.f - fem::pow2(ue(3))) * fem::cos(phi);
    ue(2) = fem::sqrt(1.f - fem::pow2(ue(3))) * fem::sin(phi);
    FEM_DO_SAFE(j, 1, 3) {
      p(n + 2, j) = pa * ue(j);
      p(n + 3, j) = -pa * ue(j);
    }
    p(n + 2, 4) = fem::sqrt(fem::pow2(pa) + fem::pow2(p(n + 2, 5)));
    p(n + 3, 4) = fem::sqrt(fem::pow2(pa) + fem::pow2(p(n + 3, 5)));
    ludbrb(n + 2, n + 3, 0.f, 0.f, dpc(1) / dpc(4), dpc(2) / dpc(4),
      dpc(3) / dpc(4));
  }
  else {
    np = 0;
    FEM_DO_SAFE(i, ic1, ic2) {
      if (k(i, 1) == 1 || k(i, 1) == 2) {
        np++;
      }
    }
    ha = p(ic1, 4) * p(ic2, 4) - p(ic1, 1) * p(ic2, 1) - p(ic1, 2) * p(ic2,
      2) - p(ic1, 3) * p(ic2, 3);
    if (np >= 3 || ha <= 1.25f * p(ic1, 5) * p(ic2, 5)) {
      goto statement_260;
    }
    hd1 = 0.5f * (fem::pow2(p(n + 2, 5)) - fem::pow2(p(ic1, 5)));
    hd2 = 0.5f * (fem::pow2(p(n + 3, 5)) - fem::pow2(p(ic2, 5)));
    hr = fem::sqrt(fem::max(0.f, (fem::pow2((ha - hd1 - hd2)) -
      fem::pow2((p(n + 2, 5) * p(n + 3, 5)))) / (fem::pow2(ha) -
      fem::pow2((p(ic1, 5) * p(ic2, 5)))))) - 1.f;
    hc = fem::pow2(p(ic1, 5)) + 2.f * ha + fem::pow2(p(ic2, 5));
    hk1 = ((fem::pow2(p(ic2, 5)) + ha) * hr + hd1 - hd2) / hc;
    hk2 = ((fem::pow2(p(ic1, 5)) + ha) * hr + hd2 - hd1) / hc;
    FEM_DO_SAFE(j, 1, 4) {
      p(n + 2, j) = (1.f + hk1) * p(ic1, j) - hk2 * p(ic2, j);
      p(n + 3, j) = (1.f + hk2) * p(ic2, j) - hk1 * p(ic1, j);
    }
  }
  FEM_DO_SAFE(j, 1, 4) {
    v(n + 1, j) = v(ic1, j);
    v(n + 2, j) = v(ic1, j);
    v(n + 3, j) = v(ic2, j);
  }
  v(n + 1, 5) = 0.f;
  v(n + 2, 5) = 0.f;
  v(n + 3, 5) = 0.f;
  n += 3;
  goto statement_300;
  //C
  //C...Else form one particle from the flavours available, if possible.
  statement_260:
  k(n + 1, 5) = n + 2;
  if (fem::iabs(k(ic1, 2)) > 100 && fem::iabs(k(ic2, 2)) > 100) {
    goto statement_320;
  }
  else if (fem::iabs(k(ic1, 2)) != 21) {
    lukfdi(cmn, k(ic1, 2), k(ic2, 2), kfldmp, k(n + 2, 2));
  }
  else {
    kfln = 1 + fem::fint((2.f + parj(2)) * rlu(cmn, 0));
    lukfdi(cmn, kfln, -kfln, kfldmp, k(n + 2, 2));
  }
  if (k(n + 2, 2) == 0) {
    goto statement_260;
  }
  p(n + 2, 5) = ulmass(cmn, k(n + 2, 2));
  //C
  //C...Find parton/particle which combines to largest extra mass.
  ir = 0;
  ha = 0.f;
  FEM_DO_SAFE(mcomb, 1, 3) {
    if (ir != 0) {
      goto statement_280;
    }
    FEM_DO_SAFE(i, fem::max(1, ip), n) {
      if (k(i, 1) <= 0 || k(i, 1) > 10 || (i >= ic1 && i <= ic2 && k(i,
          1) >= 1 && k(i, 1) <= 2)) {
        goto statement_270;
      }
      if (mcomb == 1) {
        kci = lucomp(cmn, k(i, 2));
      }
      if (mcomb == 1 && kci == 0) {
        goto statement_270;
      }
      if (mcomb == 1 && kchg(kci, 2) == 0 && i <= ns) {
        goto statement_270;
      }
      if (mcomb == 2 && fem::iabs(k(i, 2)) > 10 && fem::iabs(k(i, 2)) <= 100) {
        goto statement_270;
      }
      hcr = fem::sngl(dpc(4)) * p(i, 4) - fem::sngl(dpc(1)) * p(i,
        1) - fem::sngl(dpc(2)) * p(i, 2) - fem::sngl(dpc(3)) * p(i,
        3);
      if (hcr > ha) {
        ir = i;
        ha = hcr;
      }
      statement_270:;
    }
    statement_280:;
  }
  //C
  //C...Shuffle energy and momentum to put new particle on mass shell.
  hb = fem::pow2(pecm) + ha;
  hc = fem::pow2(p(n + 2, 5)) + ha;
  hd = fem::pow2(p(ir, 5)) + ha;
  //C******************CHANGES BY HIJING************
  hk2 = 0.0f;
  if (fem::pow2(ha) - fem::pow2((pecm * p(ir, 5))) == 0.0f || hb + hd == 0.0f) {
    goto statement_285;
  }
  //C******************
  hk2 = 0.5f * (hb * fem::sqrt((fem::pow2((hb + hc)) - 4.f * (hb +
    hd) * fem::pow2(p(n + 2, 5))) / (fem::pow2(ha) - fem::pow2((pecm * p(ir,
    5))))) - (hb + hc)) / (hb + hd);
  statement_285:
  hk1 = (0.5f * (fem::pow2(p(n + 2, 5)) - fem::pow2(pecm)) + hd * hk2) / hb;
  FEM_DO_SAFE(j, 1, 4) {
    p(n + 2, j) = (1.f + hk1) * fem::sngl(dpc(j)) - hk2 * p(ir, j);
    p(ir, j) = (1.f + hk2) * p(ir, j) - hk1 * fem::sngl(dpc(j));
    v(n + 1, j) = v(ic1, j);
    v(n + 2, j) = v(ic1, j);
  }
  v(n + 1, 5) = 0.f;
  v(n + 2, 5) = 0.f;
  n += 2;
  //C
  //C...Mark collapsed system and store daughter pointers. Iterate.
  statement_300:
  FEM_DO_SAFE(i, ic1, ic2) {
    if ((k(i, 1) == 1 || k(i, 1) == 2) && kchg(lucomp(cmn, k(i, 2)), 2) != 0) {
      k(i, 1) += 10;
      if (mstu(16) != 2) {
        k(i, 4) = nsav + 1;
        k(i, 5) = nsav + 1;
      }
      else {
        k(i, 4) = nsav + 2;
        k(i, 5) = n;
      }
    }
  }
  if (n < mstu(4) - mstu(32) - 5) {
    goto statement_140;
  }
  //C
  //C...Check flavours and invariant masses in parton systems.
  statement_320:
  np = 0;
  kfn = 0;
  kqs = 0;
  FEM_DO_SAFE(j, 1, 5) {
    dps(j) = 0e0;
  }
  FEM_DO_SAFE(i, fem::max(1, ip), n) {
    if (k(i, 1) <= 0 || k(i, 1) > 10) {
      goto statement_360;
    }
    kc = lucomp(cmn, k(i, 2));
    if (kc == 0) {
      goto statement_360;
    }
    kq = kchg(kc, 2) * fem::isign(1, k(i, 2));
    if (kq == 0) {
      goto statement_360;
    }
    np++;
    if (kq != 2) {
      kfn++;
      kqs += kq;
      mstj(93) = 1;
      dps(5) += fem::dble(ulmass(cmn, k(i, 2)));
    }
    FEM_DO_SAFE(j, 1, 4) {
      dps(j) += fem::dble(p(i, j));
    }
    //C
    //Clin-4/12/01:
    //C     np: # of partons, KFN: number of quarks and diquarks,
    //C     KC=0 for color singlet system, -1 for quarks and anti-diquarks,
    //C     1 for quarks and anti-diquarks, and 2 for gluons:
    if (k(i, 1) == 1) {
      //Clin-4/12/01     end of color singlet system.
      if (np != 1 && (kfn == 1 || kfn >= 3 || kqs != 0)) {
        luerrm(cmn, 2, "(LUPREP:) unphysical flavour combination");
      }
      //C
      //Clin-4/16/01: 'jet system' should be defined as np.ne.2:
      //C        IF(NP.NE.1.AND.DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2.LT.
      //C     &  (0.9*PARJ(32)+DPS(5))**2) CALL LUERRM(3,
      //C     &  '(LUPREP:) too small mass in jet system')
      if (np != 2 && fem::pow2(dps(4)) - fem::pow2(dps(1)) -
          fem::pow2(dps(2)) - fem::pow2(dps(3)) < fem::pow2((0.9e0 *
          fem::dble(parj(32)) + dps(5)))) {
        luerrm(cmn, 3, "(LUPREP:) too small mass in jet system");
        write(6, star), "DPS(1-5),KI1-5=", dps(1), dps(2), dps(3),
          dps(4), dps(5), "*", k(i, 1), k(i, 2), k(i, 3), k(i, 4), k(i,
          5);
      }
      //C
      np = 0;
      kfn = 0;
      kqs = 0;
      FEM_DO_SAFE(j, 1, 5) {
        dps(j) = 0e0;
      }
    }
    statement_360:;
  }
  //C
}

//C
//C*********************************************************************
//C
void
luptdi(
  common& cmn,
  int const& kfl,
  float& px,
  float& py)
{
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  //
  //C
  //C...Purpose: to generate transverse momentum according to a Gaussian.
  //C
  //C...Generate p_T and azimuthal angle, gives p_x and p_y.
  int kfla = fem::iabs(kfl);
  float pt = parj(21) * fem::sqrt(-fem::log(fem::max(1e-10f, rlu(cmn, 0))));
  if (mstj(91) == 1) {
    pt = parj(22) * pt;
  }
  if (kfla == 0 && mstj(13) <= 0) {
    pt = 0.f;
  }
  float phi = paru(2) * rlu(cmn, 0);
  px = pt * fem::cos(phi);
  py = pt * fem::sin(phi);
  //C
}

//C
//C*********************************************************************
//C
void
luzdis(
  common& cmn,
  int const& kfl1,
  int const& kfl2,
  float const& pr,
  float& z)
{
  arr_cref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  //
  int kfla = fem::int0;
  int kflb = fem::int0;
  int kflh = fem::int0;
  float fa = fem::float0;
  float fb = fem::float0;
  float fc = fem::float0;
  int mc = fem::int0;
  int ma = fem::int0;
  float zmax = fem::float0;
  int mmax = fem::int0;
  float zdiv = fem::float0;
  float fint = fem::float0;
  float zdivc = fem::float0;
  float fscb = fem::float0;
  float fpre = fem::float0;
  float fval = fem::float0;
  //C
  //C...Purpose: to generate the longitudinal splitting variable z.
  //C
  //C...Check if heavy flavour fragmentation.
  kfla = fem::iabs(kfl1);
  kflb = fem::iabs(kfl2);
  kflh = kfla;
  if (kfla >= 10) {
    kflh = fem::mod(kfla / 1000, 10);
  }
  //C
  //C...Lund symmetric scaling function: determine parameters of shape.
  if (mstj(11) == 1 || (mstj(11) == 3 && kflh <= 3)) {
    fa = parj(41);
    if (mstj(91) == 1) {
      fa = parj(43);
    }
    if (kflb >= 10) {
      fa += parj(45);
    }
    fb = parj(42) * pr;
    if (mstj(91) == 1) {
      fb = parj(44) * pr;
    }
    fc = 1.f;
    if (kfla >= 10) {
      fc = fc - parj(45);
    }
    if (kflb >= 10) {
      fc += parj(45);
    }
    mc = 1;
    if (fem::abs(fc - 1.f) > 0.01f) {
      mc = 2;
    }
    //C
    //C...Determine position of maximum. Special cases for a = 0 or a = c.
    if (fa < 0.02f) {
      ma = 1;
      zmax = 1.f;
      if (fc > fb) {
        zmax = fb / fc;
      }
    }
    else if (fem::abs(fc - fa) < 0.01f) {
      ma = 2;
      zmax = fb / (fb + fc);
    }
    else {
      ma = 3;
      zmax = 0.5f * (fb + fc - fem::sqrt(fem::pow2((fb - fc)) + 4.f *
        fa * fb)) / (fc - fa);
      if (zmax > 0.99f && fb > 100.f) {
        zmax = 1.f - fa / fb;
      }
    }
    //C
    //C...Subdivide z range if distribution very peaked near endpoint.
    mmax = 2;
    if (zmax < 0.1f) {
      mmax = 1;
      zdiv = 2.75f * zmax;
      if (mc == 1) {
        fint = 1.f - fem::log(zdiv);
      }
      else {
        zdivc = fem::pow(zdiv, (1.f - fc));
        fint = 1.f + (1.f - 1.f / zdivc) / (fc - 1.f);
      }
    }
    else if (zmax > 0.85f && fb > 1.f) {
      mmax = 3;
      fscb = fem::sqrt(4.f + fem::pow2((fc / fb)));
      zdiv = fscb - 1.f / zmax - (fc / fb) * fem::log(zmax * 0.5f * (
        fscb + fc / fb));
      if (ma >= 2) {
        zdiv += (fa / fb) * fem::log(1.f - zmax);
      }
      zdiv = fem::min(zmax, fem::max(0.f, zdiv));
      fint = 1.f + fb * (1.f - zdiv);
    }
    //C
    //C...Choice of z, preweighted for peaks at low or high z.
    statement_100:
    z = rlu(cmn, 0);
    fpre = 1.f;
    if (mmax == 1) {
      if (fint * rlu(cmn, 0) <= 1.f) {
        z = zdiv * z;
      }
      else if (mc == 1) {
        z = fem::pow(zdiv, z);
        fpre = zdiv / z;
      }
      else {
        z = 1.f / fem::pow((zdivc + z * (1.f - zdivc)), (1.f / (1.f - fc)));
        fpre = fem::pow((zdiv / z), fc);
      }
    }
    else if (mmax == 3) {
      if (fint * rlu(cmn, 0) <= 1.f) {
        z = zdiv + fem::log(z) / fb;
        fpre = fem::exp(fb * (z - zdiv));
      }
      else {
        z = zdiv + z * (1.f - zdiv);
      }
    }
    //C
    //C...Weighting according to correct formula.
    if (z <= fb / (50.f + fb) || z >= 1.f) {
      goto statement_100;
    }
    fval = fem::pow((zmax / z), fc) * fem::exp(fb * (1.f / zmax - 1.f / z));
    if (ma >= 2) {
      fval = fem::pow(((1.f - z) / (1.f - zmax)), fa) * fval;
    }
    if (fval < rlu(cmn, 0) * fpre) {
      goto statement_100;
    }
    //C
    //C...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c.
  }
  else {
    fc = parj(50 + fem::max(1, kflh));
    if (mstj(91) == 1) {
      fc = parj(59);
    }
    statement_110:
    z = rlu(cmn, 0);
    if (fc >= 0.f && fc <= 1.f) {
      if (fc > rlu(cmn, 0)) {
        z = 1.f - fem::pow(z, (1.f / 3.f));
      }
    }
    else if (fc >  - 1.f) {
      if (-4.f * fc * z * fem::pow2((1.f - z)) < rlu(cmn, 0) *
          fem::pow2((fem::pow2((1.f - z)) - fc * z))) {
        goto statement_110;
      }
    }
    else {
      if (fc > 0.f) {
        z = 1.f - fem::pow(z, (1.f / fc));
      }
      if (fc < 0.f) {
        z = fem::pow(z, (-1.f / fc));
      }
    }
  }
  //C
}

//C
//C*********************************************************************
//C
void
lustrf(
  common& cmn,
  int const& ip)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  int i = fem::int0;
  int j = fem::int0;
  arr_2d<5, 5, double> dp(fem::fill0);
  int nsav = fem::int0;
  int np = fem::int0;
  int kqsum = fem::int0;
  arr_1d<5, double> dps(fem::fill0);
  arr_1d<4, int> mju(fem::fill0);
  int kc = fem::int0;
  int kq = fem::int0;
  int ntryr = fem::int0;
  float paru12 = fem::float0;
  float paru13 = fem::float0;
  int nr = fem::int0;
  float pdrmin = fem::float0;
  int i1 = fem::int0;
  float pap = fem::float0;
  float pvp = fem::float0;
  float pdr = fem::float0;
  int ir = fem::int0;
  int nrs = fem::int0;
  int ntry = fem::int0;
  int jt = fem::int0;
  arr_1d<2, int> njs(fem::fill0);
  int js = fem::int0;
  int iu = fem::int0;
  arr_1d<3, int> iju(fem::fill0);
  arr_2d<5, 5, float> pju(fem::fill0);
  float t12 = fem::float0;
  float t13 = fem::float0;
  float t23 = fem::float0;
  float t11 = fem::float0;
  float t22 = fem::float0;
  float tsq = fem::float0;
  float t1f = fem::float0;
  float t2f = fem::float0;
  arr_1d<5, float> tju(fem::fill0);
  int ista = fem::int0;
  int ns = fem::int0;
  int is = fem::int0;
  int is1 = fem::int0;
  int is2 = fem::int0;
  double dhkc = fem::double0;
  double dhks = fem::double0;
  double dhk1 = fem::double0;
  double dhk2 = fem::double0;
  int in1 = fem::int0;
  int isav = fem::int0;
  int irankj = fem::int0;
  arr_1d<2, int> ie(fem::fill0);
  arr_1d<9, int> in(fem::fill0);
  int jq = fem::int0;
  arr_1d<3, int> kfl(fem::fill0);
  arr_1d<3, float> px(fem::fill0);
  arr_1d<3, float> py(fem::fill0);
  arr_1d<3, float> gam(fem::fill0);
  double dhc12 = fem::double0;
  double dhcx1 = fem::double0;
  double dhcx2 = fem::double0;
  double dhcxx = fem::double0;
  double dhcy1 = fem::double0;
  double dhcy2 = fem::double0;
  double dhcyx = fem::double0;
  double dhcyy = fem::double0;
  arr_1d<2, float> pr(fem::fill0);
  float z = fem::float0;
  float pxp = fem::float0;
  float pyp = fem::float0;
  arr_1d<4, double> dhg(fem::fill0);
  int in2 = fem::int0;
  arr_1d<4, double> dhm(fem::fill0);
  double dhc = fem::double0;
  double dhs1 = fem::double0;
  double dhs2 = fem::double0;
  double dhs3 = fem::double0;
  arr_1d<2, int> kfjh(fem::fill0);
  arr_1d<2, int> kfjs(fem::fill0);
  int kfls = fem::int0;
  arr_2d<4, 5, float> pjs(fem::fill0);
  int nb = fem::int0;
  float w2sum = fem::float0;
  float w2ran = fem::float0;
  arr_1d<2, int> irank(fem::fill0);
  arr_1d<3, float> pmq(fem::fill0);
  int kdump = fem::int0;
  float pr3 = fem::float0;
  float zr = fem::float0;
  int in3 = fem::int0;
  int jr = fem::int0;
  float wmin = fem::float0;
  float wrem2 = fem::float0;
  int kfl1a = fem::int0;
  int kfl2a = fem::int0;
  float pw12 = fem::float0;
  int kfldmp = fem::int0;
  double dhr1 = fem::double0;
  double dhr2 = fem::double0;
  float fd = fem::float0;
  float fa = fem::float0;
  float prev = fem::float0;
  float fb = fem::float0;
  int im = fem::int0;
  //C...Purpose: to handle the fragmentation of an arbitrary colour singlet
  //C...jet system according to the Lund string fragmentation model.
  //C
  //C...Function: four-product of two vectors.
  four(i, j) = p(i, 4) * p(j, 4) - p(i, 1) * p(j, 1) - p(i, 2) * p(j,
    2) - p(i, 3) * p(j, 3);
  dfour(i, j) = dp(i, 4) * dp(j, 4) - dp(i, 1) * dp(j, 1) - dp(i,
    2) * dp(j, 2) - dp(i, 3) * dp(j, 3);
  //C
  //C...Reset counters. Identify parton system.
  mstj(91) = 0;
  nsav = n;
  np = 0;
  kqsum = 0;
  FEM_DO_SAFE(j, 1, 5) {
    dps(j) = 0e0;
  }
  mju(1) = 0;
  mju(2) = 0;
  i = ip - 1;
  statement_110:
  i++;
  if (i > fem::min(n, mstu(4) - mstu(32))) {
    luerrm(cmn, 12, "(LUSTRF:) failed to reconstruct jet system");
    if (mstu(21) >= 1) {
      return;
    }
  }
  if (k(i, 1) != 1 && k(i, 1) != 2 && k(i, 1) != 41) {
    goto statement_110;
  }
  kc = lucomp(cmn, k(i, 2));
  if (kc == 0) {
    goto statement_110;
  }
  kq = kchg(kc, 2) * fem::isign(1, k(i, 2));
  if (kq == 0) {
    goto statement_110;
  }
  if (n + 5 * np + 11 > mstu(4) - mstu(32) - 5) {
    luerrm(cmn, 11, "(LUSTRF:) no more memory left in LUJETS");
    if (mstu(21) >= 1) {
      return;
    }
  }
  //C
  //C...Take copy of partons to be considered. Check flavour sum.
  np++;
  FEM_DO_SAFE(j, 1, 5) {
    k(n + np, j) = k(i, j);
    p(n + np, j) = p(i, j);
    dps(j) += fem::dble(p(i, j));
  }
  k(n + np, 3) = i;
  if (fem::pow2(p(n + np, 4)) < fem::pow2(p(n + np, 1)) + fem::pow2(p(n + np,
      2)) + fem::pow2(p(n + np, 3))) {
    p(n + np, 4) = fem::sqrt(fem::pow2(p(n + np, 1)) + fem::pow2(p(n + np,
      2)) + fem::pow2(p(n + np, 3)) + fem::pow2(p(n + np, 5)));
    dps(4) += fem::dble(fem::max(0.f, p(n + np, 4) - p(i, 4)));
  }
  if (kq != 2) {
    kqsum += kq;
  }
  if (k(i, 1) == 41) {
    kqsum += 2 * kq;
    if (kqsum == kq) {
      mju(1) = n + np;
    }
    if (kqsum != kq) {
      mju(2) = n + np;
    }
  }
  if (k(i, 1) == 2 || k(i, 1) == 41) {
    goto statement_110;
  }
  if (kqsum != 0) {
    luerrm(cmn, 12, "(LUSTRF:) unphysical flavour combination");
    if (mstu(21) >= 1) {
      return;
    }
  }
  //C
  //C...Boost copied system to CM frame (for better numerical precision).
  ludbrb(n + 1, n + np, 0.f, 0.f, -dps(1) / dps(4), -dps(2) / dps(4),
    -dps(3) / dps(4));
  //C
  //C...Search for very nearby partons that may be recombined.
  ntryr = 0;
  paru12 = paru(12);
  paru13 = paru(13);
  mju(3) = mju(1);
  mju(4) = mju(2);
  nr = np;
  statement_130:
  if (nr >= 3) {
    pdrmin = 2.f * paru12;
    FEM_DO_SAFE(i, n + 1, n + nr) {
      if (i == n + nr && fem::iabs(k(n + 1, 2)) != 21) {
        goto statement_140;
      }
      i1 = i + 1;
      if (i == n + nr) {
        i1 = n + 1;
      }
      if (k(i, 1) == 41 || k(i1, 1) == 41) {
        goto statement_140;
      }
      if (mju(1) != 0 && i1 < mju(1) && fem::iabs(k(i1, 2)) != 21) {
        goto statement_140;
      }
      if (mju(2) != 0 && i > mju(2) && fem::iabs(k(i, 2)) != 21) {
        goto statement_140;
      }
      pap = fem::sqrt((fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) + fem::pow2(p(i,
        3))) * (fem::pow2(p(i1, 1)) + fem::pow2(p(i1, 2)) + fem::pow2(p(i1,
        3))));
      pvp = p(i, 1) * p(i1, 1) + p(i, 2) * p(i1, 2) + p(i, 3) * p(i1, 3);
      pdr = 4.f * fem::pow2((pap - pvp)) / (fem::pow2(paru13) * pap +
        2.f * (pap - pvp));
      if (pdr < pdrmin) {
        ir = i;
        pdrmin = pdr;
      }
      statement_140:;
    }
    //C
    //C...Recombine very nearby partons to avoid machine precision problems.
    if (pdrmin < paru12 && ir == n + nr) {
      FEM_DO_SAFE(j, 1, 4) {
        p(n + 1, j) += p(n + nr, j);
      }
      p(n + 1, 5) = fem::sqrt(fem::max(0.f, fem::pow2(p(n + 1, 4)) -
        fem::pow2(p(n + 1, 1)) - fem::pow2(p(n + 1, 2)) - fem::pow2(p(n + 1,
        3))));
      nr = nr - 1;
      goto statement_130;
    }
    else if (pdrmin < paru12) {
      FEM_DO_SAFE(j, 1, 4) {
        p(ir, j) += p(ir + 1, j);
      }
      p(ir, 5) = fem::sqrt(fem::max(0.f, fem::pow2(p(ir, 4)) - fem::pow2(p(ir,
        1)) - fem::pow2(p(ir, 2)) - fem::pow2(p(ir, 3))));
      FEM_DO_SAFE(i, ir + 1, n + nr - 1) {
        k(i, 2) = k(i + 1, 2);
        FEM_DO_SAFE(j, 1, 5) {
          p(i, j) = p(i + 1, j);
        }
      }
      if (ir == n + nr - 1) {
        k(ir, 2) = k(n + nr, 2);
      }
      nr = nr - 1;
      if (mju(1) > ir) {
        mju(1) = mju(1) - 1;
      }
      if (mju(2) > ir) {
        mju(2) = mju(2) - 1;
      }
      goto statement_130;
    }
  }
  ntryr++;
  //C
  //C...Reset particle counter. Skip ahead if no junctions are present;
  //C...this is usually the case!
  nrs = fem::max(5 * nr + 11, np);
  ntry = 0;
  statement_180:
  ntry++;
  if (ntry > 100 && ntryr <= 4) {
    paru12 = 4.f * paru12;
    paru13 = 2.f * paru13;
    goto statement_130;
  }
  else if (ntry > 100) {
    luerrm(cmn, 14, "(LUSTRF:) caught in infinite loop");
    if (mstu(21) >= 1) {
      return;
    }
  }
  i = n + nrs;
  if (mju(1) == 0 && mju(2) == 0) {
    goto statement_500;
  }
  FEM_DO_SAFE(jt, 1, 2) {
    njs(jt) = 0;
    if (mju(jt) == 0) {
      goto statement_490;
    }
    js = 3 - 2 * jt;
    //C
    //C...Find and sum up momentum on three sides of junction. Check flavours.
    FEM_DO_SAFE(iu, 1, 3) {
      iju(iu) = 0;
      FEM_DO_SAFE(j, 1, 5) {
        pju(iu, j) = 0.f;
      }
    }
    iu = 0;
    FEM_DOSTEP(i1, n + 1 + (jt - 1) * (nr - 1), n + nr + (jt - 1) * (1 - nr),
      js) {
      if (k(i1, 2) != 21 && iu <= 2) {
        iu++;
        iju(iu) = i1;
      }
      FEM_DO_SAFE(j, 1, 4) {
        pju(iu, j) += p(i1, j);
      }
    }
    FEM_DO_SAFE(iu, 1, 3) {
      pju(iu, 5) = fem::sqrt(fem::pow2(pju(iu, 1)) + fem::pow2(pju(iu,
        2)) + fem::pow2(pju(iu, 3)));
    }
    if (k(iju(3), 2) / 100 != 10 * k(iju(1), 2) + k(iju(2), 2) && k(iju(3),
        2) / 100 != 10 * k(iju(2), 2) + k(iju(1), 2)) {
      luerrm(cmn, 12, "(LUSTRF:) unphysical flavour combination");
      if (mstu(21) >= 1) {
        return;
      }
    }
    //C
    //C...Calculate (approximate) boost to rest frame of junction.
    t12 = (pju(1, 1) * pju(2, 1) + pju(1, 2) * pju(2, 2) + pju(1,
      3) * pju(2, 3)) / (pju(1, 5) * pju(2, 5));
    t13 = (pju(1, 1) * pju(3, 1) + pju(1, 2) * pju(3, 2) + pju(1,
      3) * pju(3, 3)) / (pju(1, 5) * pju(3, 5));
    t23 = (pju(2, 1) * pju(3, 1) + pju(2, 2) * pju(3, 2) + pju(2,
      3) * pju(3, 3)) / (pju(2, 5) * pju(3, 5));
    t11 = fem::sqrt((2.f / 3.f) * (1.f - t12) * (1.f - t13) / (1.f - t23));
    t22 = fem::sqrt((2.f / 3.f) * (1.f - t12) * (1.f - t23) / (1.f - t13));
    tsq = fem::sqrt((2.f * t11 * t22 + t12 - 1.f) * (1.f + t12));
    t1f = (tsq - t22 * (1.f + t12)) / (1.f - fem::pow2(t12));
    t2f = (tsq - t11 * (1.f + t12)) / (1.f - fem::pow2(t12));
    FEM_DO_SAFE(j, 1, 3) {
      tju(j) = -(t1f * pju(1, j) / pju(1, 5) + t2f * pju(2, j) / pju(2, 5));
    }
    tju(4) = fem::sqrt(1.f + fem::pow2(tju(1)) + fem::pow2(tju(2)) +
      fem::pow2(tju(3)));
    FEM_DO_SAFE(iu, 1, 3) {
      pju(iu, 5) = tju(4) * pju(iu, 4) - tju(1) * pju(iu, 1) - tju(2) * pju(iu,
        2) - tju(3) * pju(iu, 3);
    }
    //C
    //C...Put junction at rest if motion could give inconsistencies.
    if (pju(1, 5) + pju(2, 5) > pju(1, 4) + pju(2, 4)) {
      FEM_DO_SAFE(j, 1, 3) {
        tju(j) = 0.f;
      }
      tju(4) = 1.f;
      pju(1, 5) = pju(1, 4);
      pju(2, 5) = pju(2, 4);
      pju(3, 5) = pju(3, 4);
    }
    //C
    //C...Start preparing for fragmentation of two strings from junction.
    ista = i;
    FEM_DO_SAFE(iu, 1, 2) {
      ns = iju(iu + 1) - iju(iu);
      //C
      //C...Junction strings: find longitudinal string directions.
      FEM_DO_SAFE(is, 1, ns) {
        is1 = iju(iu) + is - 1;
        is2 = iju(iu) + is;
        FEM_DO_SAFE(j, 1, 5) {
          dp(1, j) = fem::dble(0.5f * p(is1, j));
          if (is == 1) {
            dp(1, j) = fem::dble(p(is1, j));
          }
          dp(2, j) = fem::dble(0.5f * p(is2, j));
          if (is == ns) {
            dp(2, j) = -fem::dble(pju(iu, j));
          }
        }
        if (is == ns) {
          dp(2, 4) = fem::dble(fem::sqrt(fem::pow2(pju(iu, 1)) +
            fem::pow2(pju(iu, 2)) + fem::pow2(pju(iu, 3))));
        }
        if (is == ns) {
          dp(2, 5) = 0e0;
        }
        dp(3, 5) = dfour(1, 1);
        dp(4, 5) = dfour(2, 2);
        dhkc = dfour(1, 2);
        if (dp(3, 5) + 2e0 * dhkc + dp(4, 5) <= 0e0) {
          dp(1, 4) = fem::sqrt(fem::pow2(dp(1, 1)) + fem::pow2(dp(1,
            2)) + fem::pow2(dp(1, 3)));
          dp(2, 4) = fem::sqrt(fem::pow2(dp(2, 1)) + fem::pow2(dp(2,
            2)) + fem::pow2(dp(2, 3)));
          dp(3, 5) = 0e0;
          dp(4, 5) = 0e0;
          dhkc = dfour(1, 2);
        }
        dhks = fem::sqrt(fem::pow2(dhkc) - dp(3, 5) * dp(4, 5));
        dhk1 = 0.5e0 * ((dp(4, 5) + dhkc) / dhks - 1e0);
        dhk2 = 0.5e0 * ((dp(3, 5) + dhkc) / dhks - 1e0);
        in1 = n + nr + 4 * is - 3;
        p(in1, 5) = fem::sngl(fem::sqrt(dp(3, 5) + 2e0 * dhkc + dp(4, 5)));
        FEM_DO_SAFE(j, 1, 4) {
          p(in1, j) = fem::sngl((1e0 + dhk1) * dp(1, j) - dhk2 * dp(2, j));
          p(in1 + 1, j) = fem::sngl((1e0 + dhk2) * dp(2, j) - dhk1 * dp(1, j));
        }
      }
      //C
      //C...Junction strings: initialize flavour, momentum and starting pos.
      isav = i;
      statement_270:
      ntry++;
      if (ntry > 100 && ntryr <= 4) {
        paru12 = 4.f * paru12;
        paru13 = 2.f * paru13;
        goto statement_130;
      }
      else if (ntry > 100) {
        luerrm(cmn, 14, "(LUSTRF:) caught in infinite loop");
        if (mstu(21) >= 1) {
          return;
        }
      }
      i = isav;
      irankj = 0;
      ie(1) = k(n + 1 + (jt / 2) * (np - 1), 3);
      in(4) = n + nr + 1;
      in(5) = in(4) + 1;
      in(6) = n + nr + 4 * ns + 1;
      FEM_DO_SAFE(jq, 1, 2) {
        FEM_DOSTEP(in1, n + nr + 2 + jq, n + nr + 4 * ns - 2 + jq, 4) {
          p(in1, 1) = 2 - jq;
          p(in1, 2) = jq - 1;
          p(in1, 3) = 1.f;
        }
      }
      kfl(1) = k(iju(iu), 2);
      px(1) = 0.f;
      py(1) = 0.f;
      gam(1) = 0.f;
      FEM_DO_SAFE(j, 1, 5) {
        pju(iu + 3, j) = 0.f;
      }
      //C
      //C...Junction strings: find initial transverse directions.
      FEM_DO_SAFE(j, 1, 4) {
        dp(1, j) = fem::dble(p(in(4), j));
        dp(2, j) = fem::dble(p(in(4) + 1, j));
        dp(3, j) = 0e0;
        dp(4, j) = 0e0;
      }
      dp(1, 4) = fem::sqrt(fem::pow2(dp(1, 1)) + fem::pow2(dp(1,
        2)) + fem::pow2(dp(1, 3)));
      dp(2, 4) = fem::sqrt(fem::pow2(dp(2, 1)) + fem::pow2(dp(2,
        2)) + fem::pow2(dp(2, 3)));
      dp(5, 1) = dp(1, 1) / dp(1, 4) - dp(2, 1) / dp(2, 4);
      dp(5, 2) = dp(1, 2) / dp(1, 4) - dp(2, 2) / dp(2, 4);
      dp(5, 3) = dp(1, 3) / dp(1, 4) - dp(2, 3) / dp(2, 4);
      if (fem::pow2(dp(5, 1)) <= fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
        dp(3, 1) = 1e0;
      }
      if (fem::pow2(dp(5, 1)) > fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
        dp(3, 3) = 1e0;
      }
      if (fem::pow2(dp(5, 2)) <= fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
        dp(4, 2) = 1e0;
      }
      if (fem::pow2(dp(5, 2)) > fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
        dp(4, 3) = 1e0;
      }
      dhc12 = dfour(1, 2);
      dhcx1 = dfour(3, 1) / dhc12;
      dhcx2 = dfour(3, 2) / dhc12;
      dhcxx = 1e0 / fem::sqrt(1e0 + 2e0 * dhcx1 * dhcx2 * dhc12);
      dhcy1 = dfour(4, 1) / dhc12;
      dhcy2 = dfour(4, 2) / dhc12;
      dhcyx = dhcxx * (dhcx1 * dhcy2 + dhcx2 * dhcy1) * dhc12;
      dhcyy = 1e0 / fem::sqrt(1e0 + 2e0 * dhcy1 * dhcy2 * dhc12 -
        fem::pow2(dhcyx));
      FEM_DO_SAFE(j, 1, 4) {
        dp(3, j) = dhcxx * (dp(3, j) - dhcx2 * dp(1, j) - dhcx1 * dp(2, j));
        p(in(6), j) = fem::sngl(dp(3, j));
        p(in(6) + 1, j) = fem::sngl(dhcyy * (dp(4, j) - dhcy2 * dp(1,
          j) - dhcy1 * dp(2, j) - dhcyx * dp(3, j)));
      }
      //C
      //C...Junction strings: produce new particle, origin.
      statement_320:
      i++;
      if (2 * i - nsav >= mstu(4) - mstu(32) - 5) {
        luerrm(cmn, 11, "(LUSTRF:) no more memory left in LUJETS");
        if (mstu(21) >= 1) {
          return;
        }
      }
      irankj++;
      k(i, 1) = 1;
      k(i, 3) = ie(1);
      k(i, 4) = 0;
      k(i, 5) = 0;
      //C
      //C...Junction strings: generate flavour, hadron, pT, z and Gamma.
      statement_330:
      lukfdi(cmn, kfl(1), 0, kfl(3), k(i, 2));
      if (k(i, 2) == 0) {
        goto statement_270;
      }
      if (mstj(12) >= 3 && irankj == 1 && fem::iabs(kfl(1)) <= 10 &&
          fem::iabs(kfl(3)) > 10) {
        if (rlu(cmn, 0) > parj(19)) {
          goto statement_330;
        }
      }
      p(i, 5) = ulmass(cmn, k(i, 2));
      luptdi(cmn, kfl(1), px(3), py(3));
      pr(1) = fem::pow2(p(i, 5)) + fem::pow2((px(1) + px(3))) +
        fem::pow2((py(1) + py(3)));
      luzdis(cmn, kfl(1), kfl(3), pr(1), z);
      gam(3) = (1.f - z) * (gam(1) + pr(1) / z);
      FEM_DO_SAFE(j, 1, 3) {
        in(j) = in(3 + j);
      }
      //C
      //C...Junction strings: stepping within or from 'low' string region easy.
      if (in(1) + 1 == in(2) && z * p(in(1) + 2, 3) * p(in(2) + 2,
          3) * fem::pow2(p(in(1), 5)) >= pr(1)) {
        p(in(1) + 2, 4) = z * p(in(1) + 2, 3);
        p(in(2) + 2, 4) = pr(1) / (p(in(1) + 2, 4) * fem::pow2(p(in(1), 5)));
        FEM_DO_SAFE(j, 1, 4) {
          p(i, j) = (px(1) + px(3)) * p(in(3), j) + (py(1) + py(3)) *
            p(in(3) + 1, j);
        }
        goto statement_420;
      }
      else if (in(1) + 1 == in(2)) {
        p(in(2) + 2, 4) = p(in(2) + 2, 3);
        p(in(2) + 2, 1) = 1.f;
        in(2) += 4;
        if (in(2) > n + nr + 4 * ns) {
          goto statement_270;
        }
        if (four(in(1), in(2)) <= 1e-2f) {
          p(in(1) + 2, 4) = p(in(1) + 2, 3);
          p(in(1) + 2, 1) = 0.f;
          in(1) += 4;
        }
      }
      //C
      //C...Junction strings: find new transverse directions.
      statement_360:
      if (in(1) > n + nr + 4 * ns || in(2) > n + nr + 4 * ns || in(1) > in(2)) {
        goto statement_270;
      }
      if (in(1) != in(4) || in(2) != in(5)) {
        FEM_DO_SAFE(j, 1, 4) {
          dp(1, j) = fem::dble(p(in(1), j));
          dp(2, j) = fem::dble(p(in(2), j));
          dp(3, j) = 0e0;
          dp(4, j) = 0e0;
        }
        dp(1, 4) = fem::sqrt(fem::pow2(dp(1, 1)) + fem::pow2(dp(1,
          2)) + fem::pow2(dp(1, 3)));
        dp(2, 4) = fem::sqrt(fem::pow2(dp(2, 1)) + fem::pow2(dp(2,
          2)) + fem::pow2(dp(2, 3)));
        dhc12 = dfour(1, 2);
        //Clin-5/2012:
        //C        IF(DHC12.LE.1E-2) THEN
        if (dhc12 <= 1e-2) {
          p(in(1) + 2, 4) = p(in(1) + 2, 3);
          p(in(1) + 2, 1) = 0.f;
          in(1) += 4;
          goto statement_360;
        }
        in(3) = n + nr + 4 * ns + 5;
        dp(5, 1) = dp(1, 1) / dp(1, 4) - dp(2, 1) / dp(2, 4);
        dp(5, 2) = dp(1, 2) / dp(1, 4) - dp(2, 2) / dp(2, 4);
        dp(5, 3) = dp(1, 3) / dp(1, 4) - dp(2, 3) / dp(2, 4);
        if (fem::pow2(dp(5, 1)) <= fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
          dp(3, 1) = 1e0;
        }
        if (fem::pow2(dp(5, 1)) > fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
          dp(3, 3) = 1e0;
        }
        if (fem::pow2(dp(5, 2)) <= fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
          dp(4, 2) = 1e0;
        }
        if (fem::pow2(dp(5, 2)) > fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
          dp(4, 3) = 1e0;
        }
        dhcx1 = dfour(3, 1) / dhc12;
        dhcx2 = dfour(3, 2) / dhc12;
        dhcxx = 1e0 / fem::sqrt(1e0 + 2e0 * dhcx1 * dhcx2 * dhc12);
        dhcy1 = dfour(4, 1) / dhc12;
        dhcy2 = dfour(4, 2) / dhc12;
        dhcyx = dhcxx * (dhcx1 * dhcy2 + dhcx2 * dhcy1) * dhc12;
        dhcyy = 1e0 / fem::sqrt(1e0 + 2e0 * dhcy1 * dhcy2 * dhc12 -
          fem::pow2(dhcyx));
        FEM_DO_SAFE(j, 1, 4) {
          dp(3, j) = dhcxx * (dp(3, j) - dhcx2 * dp(1, j) - dhcx1 * dp(2, j));
          p(in(3), j) = fem::sngl(dp(3, j));
          p(in(3) + 1, j) = fem::sngl(dhcyy * (dp(4, j) - dhcy2 * dp(1,
            j) - dhcy1 * dp(2, j) - dhcyx * dp(3, j)));
        }
        //C...Express pT with respect to new axes, if sensible.
        pxp = -(px(3) * four(in(6), in(3)) + py(3) * four(in(6) + 1, in(3)));
        pyp = -(px(3) * four(in(6), in(3) + 1) + py(3) * four(in(6) + 1,
          in(3) + 1));
        if (fem::abs(fem::pow2(pxp) + fem::pow2(pyp) - fem::pow2(px(
            3)) - fem::pow2(py(3))) < 0.01f) {
          px(3) = pxp;
          py(3) = pyp;
        }
      }
      //C
      //C...Junction strings: sum up known four-momentum, coefficients for m2.
      FEM_DO_SAFE(j, 1, 4) {
        dhg(j) = 0e0;
        p(i, j) = px(1) * p(in(6), j) + py(1) * p(in(6) + 1, j) + px(
          3) * p(in(3), j) + py(3) * p(in(3) + 1, j);
        FEM_DOSTEP(in1, in(4), in(1) - 4, 4) {
          p(i, j) += p(in1 + 2, 3) * p(in1, j);
        }
        FEM_DOSTEP(in2, in(5), in(2) - 4, 4) {
          p(i, j) += p(in2 + 2, 3) * p(in2, j);
        }
      }
      dhm(1) = fem::dble(four(i, i));
      dhm(2) = fem::dble(2.f * four(i, in(1)));
      dhm(3) = fem::dble(2.f * four(i, in(2)));
      dhm(4) = fem::dble(2.f * four(in(1), in(2)));
      //C
      //C...Junction strings: find coefficients for Gamma expression.
      FEM_DOSTEP(in2, in(1) + 1, in(2), 4) {
        FEM_DOSTEP(in1, in(1), in2 - 1, 4) {
          dhc = fem::dble(2.f * four(in1, in2));
          dhg(1) += fem::dble(p(in1 + 2, 1) * p(in2 + 2, 1)) * dhc;
          if (in1 == in(1)) {
            dhg(2) = dhg(2) - fem::dble(p(in2 + 2, 1)) * dhc;
          }
          if (in2 == in(2)) {
            dhg(3) += fem::dble(p(in1 + 2, 1)) * dhc;
          }
          if (in1 == in(1) && in2 == in(2)) {
            dhg(4) = dhg(4) - dhc;
          }
        }
      }
      //C
      //C...Junction strings: solve (m2, Gamma) equation system for energies.
      dhs1 = dhm(3) * dhg(4) - dhm(4) * dhg(3);
      //Clin-5/2012:
      //C      IF(ABS(DHS1).LT.1E-4) GOTO 270
      if (fem::dabs(dhs1) < 1e-4) {
        goto statement_270;
      }
      dhs2 = dhm(4) * (fem::dble(gam(3)) - dhg(1)) - dhm(2) * dhg(
        3) - dhg(4) * (fem::pow2(fem::dble(p(i, 5))) - dhm(1)) + dhg(
        2) * dhm(3);
      dhs3 = dhm(2) * (fem::dble(gam(3)) - dhg(1)) - dhg(2) * (
        fem::pow2(fem::dble(p(i, 5))) - dhm(1));
      p(in(2) + 2, 4) = 0.5f * fem::sngl(fem::sqrt(fem::max(0e0,
        fem::pow2(dhs2) - 4e0 * dhs1 * dhs3)) / fem::abs(dhs1) -
        dhs2 / dhs1);
      if (dhm(2) + dhm(4) * fem::dble(p(in(2) + 2, 4)) <= 0e0) {
        goto statement_270;
      }
      p(in(1) + 2, 4) = (fem::pow2(p(i, 5)) - fem::sngl(dhm(1)) -
        fem::sngl(dhm(3)) * p(in(2) + 2, 4)) / (fem::sngl(dhm(2)) +
        fem::sngl(dhm(4)) * p(in(2) + 2, 4));
      //C
      //C...Junction strings: step to new region if necessary.
      if (p(in(2) + 2, 4) > p(in(2) + 2, 3)) {
        p(in(2) + 2, 4) = p(in(2) + 2, 3);
        p(in(2) + 2, 1) = 1.f;
        in(2) += 4;
        if (in(2) > n + nr + 4 * ns) {
          goto statement_270;
        }
        if (four(in(1), in(2)) <= 1e-2f) {
          p(in(1) + 2, 4) = p(in(1) + 2, 3);
          p(in(1) + 2, 1) = 0.f;
          in(1) += 4;
        }
        goto statement_360;
      }
      else if (p(in(1) + 2, 4) > p(in(1) + 2, 3)) {
        p(in(1) + 2, 4) = p(in(1) + 2, 3);
        p(in(1) + 2, 1) = 0.f;
        in(1) += js;
        goto statement_710;
      }
      //C
      //C...Junction strings: particle four-momentum, remainder, loop back.
      statement_420:
      FEM_DO_SAFE(j, 1, 4) {
        p(i, j) += p(in(1) + 2, 4) * p(in(1), j) + p(in(2) + 2, 4) * p(in(2),
          j);
        pju(iu + 3, j) += p(i, j);
      }
      if (p(i, 4) <= 0.f) {
        goto statement_270;
      }
      pju(iu + 3, 5) = tju(4) * pju(iu + 3, 4) - tju(1) * pju(iu + 3,
        1) - tju(2) * pju(iu + 3, 2) - tju(3) * pju(iu + 3, 3);
      if (pju(iu + 3, 5) < pju(iu, 5)) {
        kfl(1) = -kfl(3);
        px(1) = -px(3);
        py(1) = -py(3);
        gam(1) = gam(3);
        if (in(3) != in(6)) {
          FEM_DO_SAFE(j, 1, 4) {
            p(in(6), j) = p(in(3), j);
            p(in(6) + 1, j) = p(in(3) + 1, j);
          }
        }
        FEM_DO_SAFE(jq, 1, 2) {
          in(3 + jq) = in(jq);
          p(in(jq) + 2, 3) = p(in(jq) + 2, 3) - p(in(jq) + 2, 4);
          p(in(jq) + 2, 1) = p(in(jq) + 2, 1) - (3 - 2 * jq) * p(in(jq) + 2, 4);
        }
        goto statement_320;
      }
      //C
      //C...Junction strings: save quantities left after each string.
      if (fem::iabs(kfl(1)) > 10) {
        goto statement_270;
      }
      i = i - 1;
      kfjh(iu) = kfl(1);
      FEM_DO_SAFE(j, 1, 4) {
        pju(iu + 3, j) = pju(iu + 3, j) - p(i + 1, j);
      }
    }
    //C
    //C...Junction strings: put together to new effective string endpoint.
    njs(jt) = i - ista;
    kfjs(jt) = k(k(mju(jt + 2), 3), 2);
    kfls = 2 * fem::fint(rlu(cmn, 0) + 3.f * parj(4) / (1.f + 3.f *
      parj(4))) + 1;
    if (kfjh(1) == kfjh(2)) {
      kfls = 3;
    }
    if (ista != i) {
      kfjs(jt) = fem::isign(1000 * fem::max(fem::iabs(kfjh(1)),
        fem::iabs(kfjh(2))) + 100 * fem::min(fem::iabs(kfjh(1)),
        fem::iabs(kfjh(2))) + kfls, kfjh(1));
    }
    FEM_DO_SAFE(j, 1, 4) {
      pjs(jt, j) = pju(1, j) + pju(2, j) + p(mju(jt), j);
      pjs(jt + 2, j) = pju(4, j) + pju(5, j);
    }
    pjs(jt, 5) = fem::sqrt(fem::max(0.f, fem::pow2(pjs(jt, 4)) -
      fem::pow2(pjs(jt, 1)) - fem::pow2(pjs(jt, 2)) - fem::pow2(pjs(jt,
      3))));
    statement_490:;
  }
  //C
  //C...Open versus closed strings. Choose breakup region for latter.
  statement_500:
  if (mju(1) != 0 && mju(2) != 0) {
    ns = mju(2) - mju(1);
    nb = mju(1) - n;
  }
  else if (mju(1) != 0) {
    ns = n + nr - mju(1);
    nb = mju(1) - n;
  }
  else if (mju(2) != 0) {
    ns = mju(2) - n;
    nb = 1;
  }
  else if (fem::iabs(k(n + 1, 2)) != 21) {
    ns = nr - 1;
    nb = 1;
  }
  else {
    ns = nr + 1;
    w2sum = 0.f;
    FEM_DO_SAFE(is, 1, nr) {
      p(n + nr + is, 1) = 0.5f * four(n + is, n + is + 1 - nr * (is / nr));
      w2sum += p(n + nr + is, 1);
    }
    w2ran = rlu(cmn, 0) * w2sum;
    nb = 0;
    statement_520:
    nb++;
    w2sum = w2sum - p(n + nr + nb, 1);
    if (w2sum > w2ran && nb < nr) {
      goto statement_520;
    }
  }
  //C
  //C...Find longitudinal string directions (i.e. lightlike four-vectors).
  FEM_DO_SAFE(is, 1, ns) {
    is1 = n + is + nb - 1 - nr * ((is + nb - 2) / nr);
    is2 = n + is + nb - nr * ((is + nb - 1) / nr);
    FEM_DO_SAFE(j, 1, 5) {
      dp(1, j) = fem::dble(p(is1, j));
      if (fem::iabs(k(is1, 2)) == 21) {
        dp(1, j) = 0.5e0 * dp(1, j);
      }
      if (is1 == mju(1)) {
        dp(1, j) = fem::dble(pjs(1, j) - pjs(3, j));
      }
      dp(2, j) = fem::dble(p(is2, j));
      if (fem::iabs(k(is2, 2)) == 21) {
        dp(2, j) = 0.5e0 * dp(2, j);
      }
      if (is2 == mju(2)) {
        dp(2, j) = fem::dble(pjs(2, j) - pjs(4, j));
      }
    }
    dp(3, 5) = dfour(1, 1);
    dp(4, 5) = dfour(2, 2);
    dhkc = dfour(1, 2);
    if (dp(3, 5) + 2.e0 * dhkc + dp(4, 5) <= 0.e0) {
      dp(3, 5) = fem::pow2(dp(1, 5));
      dp(4, 5) = fem::pow2(dp(2, 5));
      dp(1, 4) = fem::sqrt(fem::pow2(dp(1, 1)) + fem::pow2(dp(1,
        2)) + fem::pow2(dp(1, 3)) + fem::pow2(dp(1, 5)));
      dp(2, 4) = fem::sqrt(fem::pow2(dp(2, 1)) + fem::pow2(dp(2,
        2)) + fem::pow2(dp(2, 3)) + fem::pow2(dp(2, 5)));
      dhkc = dfour(1, 2);
    }
    dhks = fem::sqrt(fem::pow2(dhkc) - dp(3, 5) * dp(4, 5));
    dhk1 = 0.5e0 * ((dp(4, 5) + dhkc) / dhks - 1.e0);
    dhk2 = 0.5e0 * ((dp(3, 5) + dhkc) / dhks - 1.e0);
    in1 = n + nr + 4 * is - 3;
    p(in1, 5) = fem::sqrt(fem::sngl(dp(3, 5) + 2.e0 * dhkc + dp(4, 5)));
    FEM_DO_SAFE(j, 1, 4) {
      p(in1, j) = fem::sngl((1.e0 + dhk1) * dp(1, j) - dhk2 * dp(2, j));
      p(in1 + 1, j) = fem::sngl((1.e0 + dhk2) * dp(2, j) - dhk1 * dp(1, j));
    }
  }
  //C
  //C...Begin initialization: sum up energy, set starting position.
  isav = i;
  statement_550:
  ntry++;
  if (ntry > 100 && ntryr <= 4) {
    paru12 = 4.f * paru12;
    paru13 = 2.f * paru13;
    goto statement_130;
  }
  else if (ntry > 100) {
    luerrm(cmn, 14, "(LUSTRF:) caught in infinite loop");
    if (mstu(21) >= 1) {
      return;
    }
  }
  i = isav;
  FEM_DO_SAFE(j, 1, 4) {
    p(n + nrs, j) = 0.f;
    FEM_DO_SAFE(is, 1, nr) {
      p(n + nrs, j) += p(n + is, j);
    }
  }
  FEM_DO_SAFE(jt, 1, 2) {
    irank(jt) = 0;
    if (mju(jt) != 0) {
      irank(jt) = njs(jt);
    }
    if (ns > nr) {
      irank(jt) = 1;
    }
    ie(jt) = k(n + 1 + (jt / 2) * (np - 1), 3);
    in(3 * jt + 1) = n + nr + 1 + 4 * (jt / 2) * (ns - 1);
    in(3 * jt + 2) = in(3 * jt + 1) + 1;
    in(3 * jt + 3) = n + nr + 4 * ns + 2 * jt - 1;
    FEM_DOSTEP(in1, n + nr + 2 + jt, n + nr + 4 * ns - 2 + jt, 4) {
      p(in1, 1) = 2 - jt;
      p(in1, 2) = jt - 1;
      p(in1, 3) = 1.f;
    }
  }
  //C
  //C...Initialize flavour and pT variables for open string.
  if (ns < nr) {
    px(1) = 0.f;
    py(1) = 0.f;
    if (ns == 1 && mju(1) + mju(2) == 0) {
      luptdi(cmn, 0, px(1), py(1));
    }
    px(2) = -px(1);
    py(2) = -py(1);
    FEM_DO_SAFE(jt, 1, 2) {
      kfl(jt) = k(ie(jt), 2);
      if (mju(jt) != 0) {
        kfl(jt) = kfjs(jt);
      }
      mstj(93) = 1;
      pmq(jt) = ulmass(cmn, kfl(jt));
      gam(jt) = 0.f;
    }
    //C
    //C...Closed string: random initial breakup flavour, pT and vertex.
  }
  else {
    kfl(3) = fem::fint(1.f + (2.f + parj(2)) * rlu(cmn, 0)) * fem::pow((-1),
      fem::fint(rlu(cmn, 0) + 0.5f));
    lukfdi(cmn, kfl(3), 0, kfl(1), kdump);
    kfl(2) = -kfl(1);
    if (fem::iabs(kfl(1)) > 10 && rlu(cmn, 0) > 0.5f) {
      kfl(2) = -(kfl(1) + fem::isign(10000, kfl(1)));
    }
    else if (fem::iabs(kfl(1)) > 10) {
      kfl(1) = -(kfl(2) + fem::isign(10000, kfl(2)));
    }
    luptdi(cmn, kfl(1), px(1), py(1));
    px(2) = -px(1);
    py(2) = -py(1);
    pr3 = fem::min(25.f, 0.1f * fem::pow2(p(n + nr + 1, 5)));
    statement_590:
    luzdis(cmn, kfl(1), kfl(2), pr3, z);
    zr = pr3 / (z * fem::pow2(p(n + nr + 1, 5)));
    if (zr >= 1.f) {
      goto statement_590;
    }
    //C
    FEM_DO_SAFE(jt, 1, 2) {
      mstj(93) = 1;
      pmq(jt) = ulmass(cmn, kfl(jt));
      gam(jt) = pr3 * (1.f - z) / z;
      in1 = n + nr + 3 + 4 * (jt / 2) * (ns - 1);
      p(in1, jt) = 1.f - z;
      p(in1, 3 - jt) = jt - 1;
      p(in1, 3) = (2 - jt) * (1.f - z) + (jt - 1) * z;
      p(in1 + 1, jt) = zr;
      p(in1 + 1, 3 - jt) = 2 - jt;
      p(in1 + 1, 3) = (2 - jt) * (1.f - zr) + (jt - 1) * zr;
    }
  }
  //C
  //C...Find initial transverse directions (i.e. spacelike four-vectors).
  FEM_DO_SAFE(jt, 1, 2) {
    if (jt == 1 || ns == nr - 1) {
      in1 = in(3 * jt + 1);
      in3 = in(3 * jt + 3);
      FEM_DO_SAFE(j, 1, 4) {
        dp(1, j) = fem::dble(p(in1, j));
        dp(2, j) = fem::dble(p(in1 + 1, j));
        dp(3, j) = 0.e0;
        dp(4, j) = 0.e0;
      }
      dp(1, 4) = fem::dsqrt(fem::pow2(dp(1, 1)) + fem::pow2(dp(1,
        2)) + fem::pow2(dp(1, 3)));
      dp(2, 4) = fem::dsqrt(fem::pow2(dp(2, 1)) + fem::pow2(dp(2,
        2)) + fem::pow2(dp(2, 3)));
      dp(5, 1) = dp(1, 1) / dp(1, 4) - dp(2, 1) / dp(2, 4);
      dp(5, 2) = dp(1, 2) / dp(1, 4) - dp(2, 2) / dp(2, 4);
      dp(5, 3) = dp(1, 3) / dp(1, 4) - dp(2, 3) / dp(2, 4);
      if (fem::pow2(dp(5, 1)) <= fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
        dp(3, 1) = 1.e0;
      }
      if (fem::pow2(dp(5, 1)) > fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
        dp(3, 3) = 1.e0;
      }
      if (fem::pow2(dp(5, 2)) <= fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
        dp(4, 2) = 1.e0;
      }
      if (fem::pow2(dp(5, 2)) > fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
        dp(4, 3) = 1.e0;
      }
      dhc12 = dfour(1, 2);
      dhcx1 = dfour(3, 1) / dhc12;
      dhcx2 = dfour(3, 2) / dhc12;
      dhcxx = 1e0 / fem::sqrt(1e0 + 2e0 * dhcx1 * dhcx2 * dhc12);
      dhcy1 = dfour(4, 1) / dhc12;
      dhcy2 = dfour(4, 2) / dhc12;
      dhcyx = dhcxx * (dhcx1 * dhcy2 + dhcx2 * dhcy1) * dhc12;
      dhcyy = 1e0 / fem::sqrt(1e0 + 2e0 * dhcy1 * dhcy2 * dhc12 -
        fem::pow2(dhcyx));
      FEM_DO_SAFE(j, 1, 4) {
        dp(3, j) = dhcxx * (dp(3, j) - dhcx2 * dp(1, j) - dhcx1 * dp(2, j));
        p(in3, j) = fem::sngl(dp(3, j));
        p(in3 + 1, j) = fem::sngl(dhcyy * (dp(4, j) - dhcy2 * dp(1,
          j) - dhcy1 * dp(2, j) - dhcyx * dp(3, j)));
      }
    }
    else {
      FEM_DO_SAFE(j, 1, 4) {
        p(in3 + 2, j) = p(in3, j);
        p(in3 + 3, j) = p(in3 + 1, j);
      }
    }
  }
  //C
  //C...Remove energy used up in junction string fragmentation.
  if (mju(1) + mju(2) > 0) {
    FEM_DO_SAFE(jt, 1, 2) {
      if (njs(jt) == 0) {
        goto statement_660;
      }
      FEM_DO_SAFE(j, 1, 4) {
        p(n + nrs, j) = p(n + nrs, j) - pjs(jt + 2, j);
      }
      statement_660:;
    }
  }
  //C
  //C...Produce new particle: side, origin.
  statement_670:
  i++;
  if (2 * i - nsav >= mstu(4) - mstu(32) - 5) {
    luerrm(cmn, 11, "(LUSTRF:) no more memory left in LUJETS");
    if (mstu(21) >= 1) {
      return;
    }
  }
  jt = fem::fint(1.5f + rlu(cmn, 0));
  if (fem::iabs(kfl(3 - jt)) > 10) {
    jt = 3 - jt;
  }
  jr = 3 - jt;
  js = 3 - 2 * jt;
  irank(jt)++;
  k(i, 1) = 1;
  k(i, 3) = ie(jt);
  k(i, 4) = 0;
  k(i, 5) = 0;
  //C
  //C...Generate flavour, hadron and pT.
  statement_680:
  lukfdi(cmn, kfl(jt), 0, kfl(3), k(i, 2));
  if (k(i, 2) == 0) {
    goto statement_550;
  }
  if (mstj(12) >= 3 && irank(jt) == 1 && fem::iabs(kfl(jt)) <= 10 &&
      fem::iabs(kfl(3)) > 10) {
    if (rlu(cmn, 0) > parj(19)) {
      goto statement_680;
    }
  }
  p(i, 5) = ulmass(cmn, k(i, 2));
  luptdi(cmn, kfl(jt), px(3), py(3));
  pr(jt) = fem::pow2(p(i, 5)) + fem::pow2((px(jt) + px(3))) +
    fem::pow2((py(jt) + py(3)));
  //C
  //C...Final hadrons for small invariant mass.
  mstj(93) = 1;
  pmq(3) = ulmass(cmn, kfl(3));
  wmin = parj(32 + mstj(11)) + pmq(1) + pmq(2) + parj(36) * pmq(3);
  if (fem::iabs(kfl(jt)) > 10 && fem::iabs(kfl(3)) > 10) {
    wmin = wmin - 0.5f * parj(36) * pmq(3);
  }
  wrem2 = four(n + nrs, n + nrs);
  if (wrem2 < 0.10f) {
    goto statement_550;
  }
  if (wrem2 < fem::pow2(fem::max(wmin * (1.f + (2.f * rlu(cmn, 0) -
      1.f) * parj(37)), parj(32) + pmq(1) + pmq(2)))) {
    goto statement_810;
  }
  //C
  //C...Choose z, which gives Gamma. Shift z for heavy flavours.
  luzdis(cmn, kfl(jt), kfl(3), pr(jt), z);
  //C
  kfl1a = fem::iabs(kfl(1));
  kfl2a = fem::iabs(kfl(2));
  if (fem::max(fem::mod(kfl1a, 10), fem::mod(kfl1a / 1000, 10),
      fem::mod(kfl2a, 10), fem::mod(kfl2a / 1000, 10)) >= 4) {
    pr(jr) = fem::pow2((pmq(jr) + pmq(3))) + fem::pow2((px(jr) - px(
      3))) + fem::pow2((py(jr) - py(3)));
    pw12 = fem::sqrt(fem::max(0.f, fem::pow2((wrem2 - pr(1) - pr(
      2))) - 4.f * pr(1) * pr(2)));
    z = (wrem2 + pr(jt) - pr(jr) + pw12 * (2.f * z - 1.f)) / (2.f * wrem2);
    pr(jr) = fem::pow2((pmq(jr) + parj(32 + mstj(11)))) + fem::pow2((
      px(jr) - px(3))) + fem::pow2((py(jr) - py(3)));
    if ((1.f - z) * (wrem2 - pr(jt) / z) < pr(jr)) {
      goto statement_810;
    }
  }
  gam(3) = (1.f - z) * (gam(jt) + pr(jt) / z);
  FEM_DO_SAFE(j, 1, 3) {
    in(j) = in(3 * jt + j);
  }
  //C
  //C...Stepping within or from 'low' string region easy.
  if (in(1) + 1 == in(2) && z * p(in(1) + 2, 3) * p(in(2) + 2, 3) *
      fem::pow2(p(in(1), 5)) >= pr(jt)) {
    p(in(jt) + 2, 4) = z * p(in(jt) + 2, 3);
    p(in(jr) + 2, 4) = pr(jt) / (p(in(jt) + 2, 4) * fem::pow2(p(in(1), 5)));
    FEM_DO_SAFE(j, 1, 4) {
      p(i, j) = (px(jt) + px(3)) * p(in(3), j) + (py(jt) + py(3)) * p(
        in(3) + 1, j);
    }
    goto statement_770;
  }
  else if (in(1) + 1 == in(2)) {
    p(in(jr) + 2, 4) = p(in(jr) + 2, 3);
    p(in(jr) + 2, jt) = 1.f;
    in(jr) += 4 * js;
    if (js * in(jr) > js * in(4 * jr)) {
      goto statement_550;
    }
    if (four(in(1), in(2)) <= 1e-2f) {
      p(in(jt) + 2, 4) = p(in(jt) + 2, 3);
      p(in(jt) + 2, jt) = 0.f;
      in(jt) += 4 * js;
    }
  }
  //C
  //C...Find new transverse directions (i.e. spacelike string vectors).
  statement_710:
  if (js * in(1) > js * in(3 * jr + 1) || js * in(2) > js * in(3 *
      jr + 2) || in(1) > in(2)) {
    goto statement_550;
  }
  if (in(1) != in(3 * jt + 1) || in(2) != in(3 * jt + 2)) {
    FEM_DO_SAFE(j, 1, 4) {
      dp(1, j) = fem::dble(p(in(1), j));
      dp(2, j) = fem::dble(p(in(2), j));
      dp(3, j) = 0.e0;
      dp(4, j) = 0.e0;
    }
    dp(1, 4) = fem::dsqrt(fem::pow2(dp(1, 1)) + fem::pow2(dp(1, 2)) +
      fem::pow2(dp(1, 3)));
    dp(2, 4) = fem::dsqrt(fem::pow2(dp(2, 1)) + fem::pow2(dp(2, 2)) +
      fem::pow2(dp(2, 3)));
    dhc12 = dfour(1, 2);
    //Clin-5/2012:
    //C        IF(DHC12.LE.1E-2) THEN
    if (dhc12 <= 1e-2) {
      p(in(jt) + 2, 4) = p(in(jt) + 2, 3);
      p(in(jt) + 2, jt) = 0.f;
      in(jt) += 4 * js;
      goto statement_710;
    }
    in(3) = n + nr + 4 * ns + 5;
    dp(5, 1) = dp(1, 1) / dp(1, 4) - dp(2, 1) / dp(2, 4);
    dp(5, 2) = dp(1, 2) / dp(1, 4) - dp(2, 2) / dp(2, 4);
    dp(5, 3) = dp(1, 3) / dp(1, 4) - dp(2, 3) / dp(2, 4);
    if (fem::pow2(dp(5, 1)) <= fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
      dp(3, 1) = 1.e0;
    }
    if (fem::pow2(dp(5, 1)) > fem::pow2(dp(5, 2)) + fem::pow2(dp(5, 3))) {
      dp(3, 3) = 1.e0;
    }
    if (fem::pow2(dp(5, 2)) <= fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
      dp(4, 2) = 1.e0;
    }
    if (fem::pow2(dp(5, 2)) > fem::pow2(dp(5, 1)) + fem::pow2(dp(5, 3))) {
      dp(4, 3) = 1.e0;
    }
    dhcx1 = dfour(3, 1) / dhc12;
    dhcx2 = dfour(3, 2) / dhc12;
    dhcxx = 1e0 / fem::sqrt(1e0 + 2e0 * dhcx1 * dhcx2 * dhc12);
    dhcy1 = dfour(4, 1) / dhc12;
    dhcy2 = dfour(4, 2) / dhc12;
    dhcyx = dhcxx * (dhcx1 * dhcy2 + dhcx2 * dhcy1) * dhc12;
    dhcyy = 1e0 / fem::sqrt(1e0 + 2e0 * dhcy1 * dhcy2 * dhc12 -
      fem::pow2(dhcyx));
    FEM_DO_SAFE(j, 1, 4) {
      dp(3, j) = dhcxx * (dp(3, j) - dhcx2 * dp(1, j) - dhcx1 * dp(2, j));
      p(in(3), j) = fem::sngl(dp(3, j));
      p(in(3) + 1, j) = fem::sngl(dhcyy * (dp(4, j) - dhcy2 * dp(1,
        j) - dhcy1 * dp(2, j) - dhcyx * dp(3, j)));
    }
    //C...Express pT with respect to new axes, if sensible.
    pxp = -(px(3) * four(in(3 * jt + 3), in(3)) + py(3) * four(in(3 *
      jt + 3) + 1, in(3)));
    pyp = -(px(3) * four(in(3 * jt + 3), in(3) + 1) + py(3) * four(in(
      3 * jt + 3) + 1, in(3) + 1));
    if (fem::abs(fem::pow2(pxp) + fem::pow2(pyp) - fem::pow2(px(3)) -
        fem::pow2(py(3))) < 0.01f) {
      px(3) = pxp;
      py(3) = pyp;
    }
  }
  //C
  //C...Sum up known four-momentum. Gives coefficients for m2 expression.
  FEM_DO_SAFE(j, 1, 4) {
    dhg(j) = 0.e0;
    p(i, j) = px(jt) * p(in(3 * jt + 3), j) + py(jt) * p(in(3 * jt + 3) + 1,
      j) + px(3) * p(in(3), j) + py(3) * p(in(3) + 1, j);
    FEM_DOSTEP(in1, in(3 * jt + 1), in(1) - 4 * js, 4 * js) {
      p(i, j) += p(in1 + 2, 3) * p(in1, j);
    }
    FEM_DOSTEP(in2, in(3 * jt + 2), in(2) - 4 * js, 4 * js) {
      p(i, j) += p(in2 + 2, 3) * p(in2, j);
    }
  }
  dhm(1) = fem::dble(four(i, i));
  dhm(2) = fem::dble(2.f * four(i, in(1)));
  dhm(3) = fem::dble(2.f * four(i, in(2)));
  dhm(4) = fem::dble(2.f * four(in(1), in(2)));
  //C
  //C...Find coefficients for Gamma expression.
  FEM_DOSTEP(in2, in(1) + 1, in(2), 4) {
    FEM_DOSTEP(in1, in(1), in2 - 1, 4) {
      dhc = fem::dble(2.f * four(in1, in2));
      dhg(1) += fem::dble(p(in1 + 2, jt) * p(in2 + 2, jt)) * dhc;
      if (in1 == in(1)) {
        dhg(2) = dhg(2) - fem::dble(fem::ffloat(js) * p(in2 + 2, jt)) * dhc;
      }
      if (in2 == in(2)) {
        dhg(3) += fem::dble(fem::ffloat(js) * p(in1 + 2, jt)) * dhc;
      }
      if (in1 == in(1) && in2 == in(2)) {
        dhg(4) = dhg(4) - dhc;
      }
    }
  }
  //C
  //C...Solve (m2, Gamma) equation system for energies taken.
  dhs1 = dhm(jr + 1) * dhg(4) - dhm(4) * dhg(jr + 1);
  //Clin-5/2012:
  //C      IF(ABS(DHS1).LT.1E-4) GOTO 550
  if (fem::dabs(dhs1) < 1e-4) {
    goto statement_550;
  }
  dhs2 = dhm(4) * (fem::dble(gam(3)) - dhg(1)) - dhm(jt + 1) * dhg(
    jr + 1) - dhg(4) * (fem::pow2(fem::dble(p(i, 5))) - dhm(1)) + dhg(
    jt + 1) * dhm(jr + 1);
  dhs3 = dhm(jt + 1) * (fem::dble(gam(3)) - dhg(1)) - dhg(jt + 1) * (
    fem::pow2(fem::dble(p(i, 5))) - dhm(1));
  p(in(jr) + 2, 4) = 0.5f * fem::sngl((fem::sqrt(fem::max(0e0,
    fem::pow2(dhs2) - 4.e0 * dhs1 * dhs3))) / fem::abs(dhs1) - dhs2 /
    dhs1);
  if (dhm(jt + 1) + dhm(4) * fem::dble(p(in(jr) + 2, 4)) <= 0.e0) {
    goto statement_550;
  }
  p(in(jt) + 2, 4) = (fem::pow2(p(i, 5)) - fem::sngl(dhm(1)) -
    fem::sngl(dhm(jr + 1)) * p(in(jr) + 2, 4)) / (fem::sngl(dhm(jt +
    1)) + fem::sngl(dhm(4)) * p(in(jr) + 2, 4));
  //C
  //C...Step to new region if necessary.
  if (p(in(jr) + 2, 4) > p(in(jr) + 2, 3)) {
    p(in(jr) + 2, 4) = p(in(jr) + 2, 3);
    p(in(jr) + 2, jt) = 1.f;
    in(jr) += 4 * js;
    if (js * in(jr) > js * in(4 * jr)) {
      goto statement_550;
    }
    if (four(in(1), in(2)) <= 1e-2f) {
      p(in(jt) + 2, 4) = p(in(jt) + 2, 3);
      p(in(jt) + 2, jt) = 0.f;
      in(jt) += 4 * js;
    }
    goto statement_710;
  }
  else if (p(in(jt) + 2, 4) > p(in(jt) + 2, 3)) {
    p(in(jt) + 2, 4) = p(in(jt) + 2, 3);
    p(in(jt) + 2, jt) = 0.f;
    in(jt) += 4 * js;
    goto statement_710;
  }
  //C
  //C...Four-momentum of particle. Remaining quantities. Loop back.
  statement_770:
  FEM_DO_SAFE(j, 1, 4) {
    p(i, j) += p(in(1) + 2, 4) * p(in(1), j) + p(in(2) + 2, 4) * p(in(2), j);
    p(n + nrs, j) = p(n + nrs, j) - p(i, j);
  }
  if (p(i, 4) <= 0.f) {
    goto statement_550;
  }
  kfl(jt) = -kfl(3);
  pmq(jt) = pmq(3);
  px(jt) = -px(3);
  py(jt) = -py(3);
  gam(jt) = gam(3);
  if (in(3) != in(3 * jt + 3)) {
    FEM_DO_SAFE(j, 1, 4) {
      p(in(3 * jt + 3), j) = p(in(3), j);
      p(in(3 * jt + 3) + 1, j) = p(in(3) + 1, j);
    }
  }
  FEM_DO_SAFE(jq, 1, 2) {
    in(3 * jt + jq) = in(jq);
    p(in(jq) + 2, 3) = p(in(jq) + 2, 3) - p(in(jq) + 2, 4);
    p(in(jq) + 2, jt) = p(in(jq) + 2, jt) - js * (3 - 2 * jq) * p(in(jq) + 2,
      4);
  }
  goto statement_670;
  //C
  //C...Final hadron: side, flavour, hadron, mass.
  statement_810:
  i++;
  k(i, 1) = 1;
  k(i, 3) = ie(jr);
  k(i, 4) = 0;
  k(i, 5) = 0;
  lukfdi(cmn, kfl(jr), -kfl(3), kfldmp, k(i, 2));
  if (k(i, 2) == 0) {
    goto statement_550;
  }
  p(i, 5) = ulmass(cmn, k(i, 2));
  pr(jr) = fem::pow2(p(i, 5)) + fem::pow2((px(jr) - px(3))) +
    fem::pow2((py(jr) - py(3)));
  //C
  //C...Final two hadrons: find common setup of four-vectors.
  jq = 1;
  if (p(in(4) + 2, 3) * p(in(5) + 2, 3) * four(in(4), in(5)) < p(in(7),
      3) * p(in(8), 3) * four(in(7), in(8))) {
    jq = 2;
  }
  dhc12 = fem::dble(four(in(3 * jq + 1), in(3 * jq + 2)));
  dhr1 = fem::dble(four(n + nrs, in(3 * jq + 2))) / dhc12;
  dhr2 = fem::dble(four(n + nrs, in(3 * jq + 1))) / dhc12;
  if (in(4) != in(7) || in(5) != in(8)) {
    px(3 - jq) = -four(n + nrs, in(3 * jq + 3)) - px(jq);
    py(3 - jq) = -four(n + nrs, in(3 * jq + 3) + 1) - py(jq);
    pr(3 - jq) = fem::pow2(p(i + fem::pow2((jt + jq - 3)) - 1, 5)) +
      fem::pow2((px(3 - jq) + (2 * jq - 3) * js * px(3))) + fem::pow2(
      (py(3 - jq) + (2 * jq - 3) * js * py(3)));
  }
  //C
  //C...Solve kinematics for final two hadrons, if possible.
  wrem2 += fem::pow2((px(1) + px(2))) + fem::pow2((py(1) + py(2)));
  fd = (fem::sqrt(pr(1)) + fem::sqrt(pr(2))) / fem::sqrt(wrem2);
  if (mju(1) + mju(2) != 0 && i == isav + 2 && fd >= 1.f) {
    goto statement_180;
  }
  if (fd >= 1.f) {
    goto statement_550;
  }
  fa = wrem2 + pr(jt) - pr(jr);
  if (mstj(11) == 2) {
    prev = 0.5f * fem::pow(fd, parj(37 + mstj(11)));
  }
  if (mstj(11) != 2) {
    prev = 0.5f * fem::exp(fem::max(-100.f, fem::log(fd) * parj(37 +
      mstj(11)) * fem::pow2((pr(1) + pr(2)))));
  }
  fb = fem::sign(fem::sqrt(fem::max(0.f, fem::pow2(fa) - 4.f *
    wrem2 * pr(jt))), js * (rlu(cmn, 0) - prev));
  kfl1a = fem::iabs(kfl(1));
  kfl2a = fem::iabs(kfl(2));
  if (fem::max(fem::mod(kfl1a, 10), fem::mod(kfl1a / 1000, 10),
      fem::mod(kfl2a, 10), fem::mod(kfl2a / 1000, 10)) >= 6) {
    fb = fem::sign(fem::sqrt(fem::max(0.f, fem::pow2(fa) - 4.f *
      wrem2 * pr(jt))), fem::ffloat(js));
  }
  FEM_DO_SAFE(j, 1, 4) {
    p(i - 1, j) = (px(jt) + px(3)) * p(in(3 * jq + 3), j) + (py(jt) +
      py(3)) * p(in(3 * jq + 3) + 1, j) + 0.5f * (fem::sngl(dhr1) * (
      fa + fb) * p(in(3 * jq + 1), j) + fem::sngl(dhr2) * (fa - fb) *
      p(in(3 * jq + 2), j)) / wrem2;
    p(i, j) = p(n + nrs, j) - p(i - 1, j);
  }
  //C
  //C...Mark jets as fragmented and give daughter pointers.
  n = i - nrs + 1;
  FEM_DO_SAFE(i, nsav + 1, nsav + np) {
    im = k(i, 3);
    k(im, 1) += 10;
    if (mstu(16) != 2) {
      k(im, 4) = nsav + 1;
      k(im, 5) = nsav + 1;
    }
    else {
      k(im, 4) = nsav + 2;
      k(im, 5) = n;
    }
  }
  //C
  //C...Document string system. Move up particles.
  nsav++;
  k(nsav, 1) = 11;
  k(nsav, 2) = 92;
  k(nsav, 3) = ip;
  k(nsav, 4) = nsav + 1;
  k(nsav, 5) = n;
  FEM_DO_SAFE(j, 1, 4) {
    p(nsav, j) = fem::sngl(dps(j));
    v(nsav, j) = v(ip, j);
  }
  p(nsav, 5) = fem::sqrt(fem::sngl(fem::max(0e0, fem::pow2(dps(4)) -
    fem::pow2(dps(1)) - fem::pow2(dps(2)) - fem::pow2(dps(3)))));
  v(nsav, 5) = 0.f;
  FEM_DO_SAFE(i, nsav + 1, n) {
    //C
    FEM_DO_SAFE(j, 1, 5) {
      k(i, j) = k(i + nrs - 1, j);
      p(i, j) = p(i + nrs - 1, j);
      v(i, j) = 0.f;
    }
  }
  //C
  //C...Order particles in rank along the chain. Update mother pointer.
  FEM_DO_SAFE(i, nsav + 1, n) {
    FEM_DO_SAFE(j, 1, 5) {
      k(i - nsav + n, j) = k(i, j);
      p(i - nsav + n, j) = p(i, j);
    }
  }
  i1 = nsav;
  FEM_DO_SAFE(i, n + 1, 2 * n - nsav) {
    if (k(i, 3) != ie(1)) {
      goto statement_880;
    }
    i1++;
    FEM_DO_SAFE(j, 1, 5) {
      k(i1, j) = k(i, j);
      p(i1, j) = p(i, j);
    }
    if (mstu(16) != 2) {
      k(i1, 3) = nsav;
    }
    statement_880:;
  }
  FEM_DOSTEP(i, 2 * n - nsav, n + 1, -1) {
    if (k(i, 3) == ie(1)) {
      goto statement_900;
    }
    i1++;
    FEM_DO_SAFE(j, 1, 5) {
      k(i1, j) = k(i, j);
      p(i1, j) = p(i, j);
    }
    if (mstu(16) != 2) {
      k(i1, 3) = nsav;
    }
    statement_900:;
  }
  //C
  //C...Boost back particle system. Set production vertices.
  ludbrb(nsav + 1, n, 0.f, 0.f, dps(1) / dps(4), dps(2) / dps(4), dps(
    3) / dps(4));
  FEM_DO_SAFE(i, nsav + 1, n) {
    //C
    FEM_DO_SAFE(j, 1, 4) {
      v(i, j) = v(ip, j);
    }
  }
  //C
}

//C
//C*********************************************************************
//C
void
luindf(
  common& cmn,
  int const& ip)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  int nsav = fem::int0;
  int njet = fem::int0;
  int kqsum = fem::int0;
  int j = fem::int0;
  arr_1d<5, double> dps(fem::fill0);
  int i = fem::int0;
  int kc = fem::int0;
  int kq = fem::int0;
  float pecm = fem::float0;
  arr_1d<3, int> nfi(fem::fill0);
  int kfa = fem::int0;
  int kfla = fem::int0;
  int kflb = fem::int0;
  int ntry = fem::int0;
  arr_1d<3, int> nfl(fem::fill0);
  arr_1d<3, int> ifet(fem::fill0);
  arr_1d<3, int> kflf(fem::fill0);
  int ip1 = fem::int0;
  int nsav1 = fem::int0;
  int kflh = fem::int0;
  arr_1d<2, int> kflo(fem::fill0);
  float wf = fem::float0;
  int nstr = fem::int0;
  arr_1d<2, float> pxo(fem::fill0);
  arr_1d<2, float> pyo(fem::fill0);
  arr_1d<2, float> wo(fem::fill0);
  int istr = fem::int0;
  int irank = fem::int0;
  int kfl1 = fem::int0;
  float px1 = fem::float0;
  float py1 = fem::float0;
  float w = fem::float0;
  int kfl2 = fem::int0;
  float px2 = fem::float0;
  float py2 = fem::float0;
  float pr = fem::float0;
  float z = fem::float0;
  float the = fem::float0;
  float phi = fem::float0;
  int kflc = fem::int0;
  int nreq = fem::int0;
  int nrem = fem::int0;
  int irem = fem::int0;
  float p2min = fem::float0;
  float p2 = fem::float0;
  int isgn = fem::int0;
  int nfet = fem::int0;
  int kflfc = fem::int0;
  int kfldmp = fem::int0;
  int kf = fem::int0;
  int npos = fem::int0;
  arr_1d<4, float> psi(fem::fill0);
  float pws = fem::float0;
  float pw = fem::float0;
  int ir1 = fem::int0;
  int ir2 = fem::int0;
  float pls = fem::float0;
  float pss = fem::float0;
  float pms = fem::float0;
  float pes = fem::float0;
  float pqs = fem::float0;
  int neco = fem::int0;
  float pfac = fem::float0;
  int i1 = fem::int0;
  //C
  //C...Purpose: to handle the fragmentation of a jet system (or a single
  //C...jet) according to independent fragmentation models.
  //C
  //C...Reset counters. Identify parton system and take copy. Check flavour.
  nsav = n;
  njet = 0;
  kqsum = 0;
  FEM_DO_SAFE(j, 1, 5) {
    dps(j) = 0.e0;
  }
  i = ip - 1;
  statement_110:
  i++;
  if (i > fem::min(n, mstu(4) - mstu(32))) {
    luerrm(cmn, 12, "(LUINDF:) failed to reconstruct jet system");
    if (mstu(21) >= 1) {
      return;
    }
  }
  if (k(i, 1) != 1 && k(i, 1) != 2) {
    goto statement_110;
  }
  kc = lucomp(cmn, k(i, 2));
  if (kc == 0) {
    goto statement_110;
  }
  kq = kchg(kc, 2) * fem::isign(1, k(i, 2));
  if (kq == 0) {
    goto statement_110;
  }
  njet++;
  if (kq != 2) {
    kqsum += kq;
  }
  FEM_DO_SAFE(j, 1, 5) {
    k(nsav + njet, j) = k(i, j);
    p(nsav + njet, j) = p(i, j);
    dps(j) += fem::dble(p(i, j));
  }
  k(nsav + njet, 3) = i;
  if (k(i, 1) == 2 || (mstj(3) <= 5 && n > i && k(i + 1, 1) == 2)) {
    goto statement_110;
  }
  if (njet != 1 && kqsum != 0) {
    luerrm(cmn, 12, "(LUINDF:) unphysical flavour combination");
    if (mstu(21) >= 1) {
      return;
    }
  }
  //C
  //C...Boost copied system to CM frame. Find CM energy and sum flavours.
  if (njet != 1) {
    ludbrb(nsav + 1, nsav + njet, 0.f, 0.f, -dps(1) / dps(4), -dps(2) / dps(4),
      -dps(3) / dps(4));
  }
  pecm = 0.f;
  FEM_DO_SAFE(j, 1, 3) {
    nfi(j) = 0;
  }
  FEM_DO_SAFE(i, nsav + 1, nsav + njet) {
    pecm += p(i, 4);
    kfa = fem::iabs(k(i, 2));
    if (kfa <= 3) {
      nfi(kfa) += fem::isign(1, k(i, 2));
    }
    else if (kfa > 1000) {
      kfla = fem::mod(kfa / 1000, 10);
      kflb = fem::mod(kfa / 100, 10);
      if (kfla <= 3) {
        nfi(kfla) += fem::isign(1, k(i, 2));
      }
      if (kflb <= 3) {
        nfi(kflb) += fem::isign(1, k(i, 2));
      }
    }
  }
  //C
  //C...Loop over attempts made. Reset counters.
  ntry = 0;
  statement_150:
  ntry++;
  n = nsav + njet;
  if (ntry > 200) {
    luerrm(cmn, 14, "(LUINDF:) caught in infinite loop");
    if (mstu(21) >= 1) {
      return;
    }
  }
  FEM_DO_SAFE(j, 1, 3) {
    nfl(j) = nfi(j);
    ifet(j) = 0;
    kflf(j) = 0;
  }
  //C
  //C...Loop over jets to be fragmented.
  FEM_DO_SAFE(ip1, nsav + 1, nsav + njet) {
    mstj(91) = 0;
    nsav1 = n;
    //C
    //C...Initial flavour and momentum values. Jet along +z axis.
    kflh = fem::iabs(k(ip1, 2));
    if (kflh > 10) {
      kflh = fem::mod(kflh / 1000, 10);
    }
    kflo(2) = 0;
    wf = p(ip1, 4) + fem::sqrt(fem::pow2(p(ip1, 1)) + fem::pow2(p(ip1,
      2)) + fem::pow2(p(ip1, 3)));
    //C
    //C...Initial values for quark or diquark jet.
    statement_170:
    if (fem::iabs(k(ip1, 2)) != 21) {
      nstr = 1;
      kflo(1) = k(ip1, 2);
      luptdi(cmn, 0, pxo(1), pyo(1));
      wo(1) = wf;
      //C
      //C...Initial values for gluon treated like random quark jet.
    }
    else if (mstj(2) <= 2) {
      nstr = 1;
      if (mstj(2) == 2) {
        mstj(91) = 1;
      }
      kflo(1) = fem::fint(1.f + (2.f + parj(2)) * rlu(cmn, 0)) * fem::pow((-1),
        fem::fint(rlu(cmn, 0) + 0.5f));
      luptdi(cmn, 0, pxo(1), pyo(1));
      wo(1) = wf;
      //C
      //C...Initial values for gluon treated like quark-antiquark jet pair,
      //C...sharing energy according to Altarelli-Parisi splitting function.
    }
    else {
      nstr = 2;
      if (mstj(2) == 4) {
        mstj(91) = 1;
      }
      kflo(1) = fem::fint(1.f + (2.f + parj(2)) * rlu(cmn, 0)) * fem::pow((-1),
        fem::fint(rlu(cmn, 0) + 0.5f));
      kflo(2) = -kflo(1);
      luptdi(cmn, 0, pxo(1), pyo(1));
      pxo(2) = -pxo(1);
      pyo(2) = -pyo(1);
      wo(1) = wf * fem::pow(rlu(cmn, 0), (1.f / 3.f));
      wo(2) = wf - wo(1);
    }
    //C
    //C...Initial values for rank, flavour, pT and W+.
    FEM_DO_SAFE(istr, 1, nstr) {
      statement_180:
      i = n;
      irank = 0;
      kfl1 = kflo(istr);
      px1 = pxo(istr);
      py1 = pyo(istr);
      w = wo(istr);
      //C
      //C...New hadron. Generate flavour and hadron species.
      statement_190:
      i++;
      if (i >= mstu(4) - mstu(32) - njet - 5) {
        luerrm(cmn, 11, "(LUINDF:) no more memory left in LUJETS");
        if (mstu(21) >= 1) {
          return;
        }
      }
      irank++;
      k(i, 1) = 1;
      k(i, 3) = ip1;
      k(i, 4) = 0;
      k(i, 5) = 0;
      statement_200:
      lukfdi(cmn, kfl1, 0, kfl2, k(i, 2));
      if (k(i, 2) == 0) {
        goto statement_180;
      }
      if (mstj(12) >= 3 && irank == 1 && fem::iabs(kfl1) <= 10 &&
          fem::iabs(kfl2) > 10) {
        if (rlu(cmn, 0) > parj(19)) {
          goto statement_200;
        }
      }
      //C
      //C...Find hadron mass. Generate four-momentum.
      p(i, 5) = ulmass(cmn, k(i, 2));
      luptdi(cmn, kfl1, px2, py2);
      p(i, 1) = px1 + px2;
      p(i, 2) = py1 + py2;
      pr = fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) + fem::pow2(p(i, 2));
      luzdis(cmn, kfl1, kfl2, pr, z);
      p(i, 3) = 0.5f * (z * w - pr / (z * w));
      p(i, 4) = 0.5f * (z * w + pr / (z * w));
      if (mstj(3) >= 1 && irank == 1 && kflh >= 4 && p(i, 3) <= 0.001f) {
        if (w >= p(i, 5) + 0.5f * parj(32)) {
          goto statement_180;
        }
        p(i, 3) = 0.0001f;
        p(i, 4) = fem::sqrt(pr);
        z = p(i, 4) / w;
      }
      //C
      //C...Remaining flavour and momentum.
      kfl1 = -kfl2;
      px1 = -px2;
      py1 = -py2;
      w = (1.f - z) * w;
      FEM_DO_SAFE(j, 1, 5) {
        v(i, j) = 0.f;
      }
      //C
      //C...Check if pL acceptable. Go back for new hadron if enough energy.
      if (mstj(3) >= 0 && p(i, 3) < 0.f) {
        i = i - 1;
      }
      if (w > parj(31)) {
        goto statement_190;
      }
      n = i;
    }
    if (fem::mod(mstj(3), 5) == 4 && n == nsav1) {
      wf += 0.1f * parj(32);
    }
    if (fem::mod(mstj(3), 5) == 4 && n == nsav1) {
      goto statement_170;
    }
    //C
    //C...Rotate jet to new direction.
    the = ulangl(cmn, p(ip1, 3), fem::sqrt(fem::pow2(p(ip1, 1)) +
      fem::pow2(p(ip1, 2))));
    phi = ulangl(cmn, p(ip1, 1), p(ip1, 2));
    ludbrb(nsav1 + 1, n, the, phi, 0e0, 0e0, 0e0);
    k(k(ip1, 3), 4) = nsav1 + 1;
    k(k(ip1, 3), 5) = n;
    //C
    //C...End of jet generation loop. Skip conservation in some cases.
  }
  if (njet == 1 || mstj(3) <= 0) {
    goto statement_470;
  }
  if (fem::mod(mstj(3), 5) != 0 && n - nsav - njet < 2) {
    goto statement_150;
  }
  //C
  //C...Subtract off produced hadron flavours, finished if zero.
  FEM_DO_SAFE(i, nsav + njet + 1, n) {
    kfa = fem::iabs(k(i, 2));
    kfla = fem::mod(kfa / 1000, 10);
    kflb = fem::mod(kfa / 100, 10);
    kflc = fem::mod(kfa / 10, 10);
    if (kfla == 0) {
      if (kflb <= 3) {
        nfl(kflb) = nfl(kflb) - fem::isign(1, k(i, 2)) * fem::pow((-1), kflb);
      }
      if (kflc <= 3) {
        nfl(kflc) += fem::isign(1, k(i, 2)) * fem::pow((-1), kflb);
      }
    }
    else {
      if (kfla <= 3) {
        nfl(kfla) = nfl(kfla) - fem::isign(1, k(i, 2));
      }
      if (kflb <= 3) {
        nfl(kflb) = nfl(kflb) - fem::isign(1, k(i, 2));
      }
      if (kflc <= 3) {
        nfl(kflc) = nfl(kflc) - fem::isign(1, k(i, 2));
      }
    }
  }
  nreq = (fem::iabs(nfl(1)) + fem::iabs(nfl(2)) + fem::iabs(nfl(3)) -
    fem::iabs(nfl(1) + nfl(2) + nfl(3))) / 2 + fem::iabs(nfl(1) + nfl(2) +
    nfl(3)) / 3;
  if (nreq == 0) {
    goto statement_320;
  }
  //C
  //C...Take away flavour of low-momentum particles until enough freedom.
  nrem = 0;
  statement_250:
  irem = 0;
  p2min = fem::pow2(pecm);
  FEM_DO_SAFE(i, nsav + njet + 1, n) {
    p2 = fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) + fem::pow2(p(i, 3));
    if (k(i, 1) == 1 && p2 < p2min) {
      irem = i;
    }
    if (k(i, 1) == 1 && p2 < p2min) {
      p2min = p2;
    }
  }
  if (irem == 0) {
    goto statement_150;
  }
  k(irem, 1) = 7;
  kfa = fem::iabs(k(irem, 2));
  kfla = fem::mod(kfa / 1000, 10);
  kflb = fem::mod(kfa / 100, 10);
  kflc = fem::mod(kfa / 10, 10);
  if (kfla >= 4 || kflb >= 4) {
    k(irem, 1) = 8;
  }
  if (k(irem, 1) == 8) {
    goto statement_250;
  }
  if (kfla == 0) {
    isgn = fem::isign(1, k(irem, 2)) * fem::pow((-1), kflb);
    if (kflb <= 3) {
      nfl(kflb) += isgn;
    }
    if (kflc <= 3) {
      nfl(kflc) = nfl(kflc) - isgn;
    }
  }
  else {
    if (kfla <= 3) {
      nfl(kfla) += fem::isign(1, k(irem, 2));
    }
    if (kflb <= 3) {
      nfl(kflb) += fem::isign(1, k(irem, 2));
    }
    if (kflc <= 3) {
      nfl(kflc) += fem::isign(1, k(irem, 2));
    }
  }
  nrem++;
  nreq = (fem::iabs(nfl(1)) + fem::iabs(nfl(2)) + fem::iabs(nfl(3)) -
    fem::iabs(nfl(1) + nfl(2) + nfl(3))) / 2 + fem::iabs(nfl(1) + nfl(2) +
    nfl(3)) / 3;
  if (nreq > nrem) {
    goto statement_250;
  }
  FEM_DO_SAFE(i, nsav + njet + 1, n) {
    if (k(i, 1) == 8) {
      k(i, 1) = 1;
    }
  }
  //C
  //C...Find combination of existing and new flavours for hadron.
  statement_280:
  nfet = 2;
  if (nfl(1) + nfl(2) + nfl(3) != 0) {
    nfet = 3;
  }
  if (nreq < nrem) {
    nfet = 1;
  }
  if (fem::iabs(nfl(1)) + fem::iabs(nfl(2)) + fem::iabs(nfl(3)) == 0) {
    nfet = 0;
  }
  FEM_DO_SAFE(j, 1, nfet) {
    ifet(j) = 1 + fem::fint((fem::iabs(nfl(1)) + fem::iabs(nfl(2)) +
      fem::iabs(nfl(3))) * rlu(cmn, 0));
    kflf(j) = fem::isign(1, nfl(1));
    if (ifet(j) > fem::iabs(nfl(1))) {
      kflf(j) = fem::isign(2, nfl(2));
    }
    if (ifet(j) > fem::iabs(nfl(1)) + fem::iabs(nfl(2))) {
      kflf(j) = fem::isign(3, nfl(3));
    }
  }
  if (nfet == 2 && (ifet(1) == ifet(2) || kflf(1) * kflf(2) > 0)) {
    goto statement_280;
  }
  if (nfet == 3 && (ifet(1) == ifet(2) || ifet(1) == ifet(3) || ifet(
      2) == ifet(3) || kflf(1) * kflf(2) < 0 || kflf(1) * kflf(
      3) < 0 || kflf(1) * (nfl(1) + nfl(2) + nfl(3)) < 0)) {
    goto statement_280;
  }
  if (nfet == 0) {
    kflf(1) = 1 + fem::fint((2.f + parj(2)) * rlu(cmn, 0));
  }
  if (nfet == 0) {
    kflf(2) = -kflf(1);
  }
  if (nfet == 1) {
    kflf(2) = fem::isign(1 + fem::fint((2.f + parj(2)) * rlu(cmn,
      0)), -kflf(1));
  }
  if (nfet <= 2) {
    kflf(3) = 0;
  }
  if (kflf(3) != 0) {
    kflfc = fem::isign(1000 * fem::max(fem::iabs(kflf(1)), fem::iabs(
      kflf(3))) + 100 * fem::min(fem::iabs(kflf(1)), fem::iabs(kflf(3))) + 1,
      kflf(1));
    if (kflf(1) == kflf(3) || (1.f + 3.f * parj(4)) * rlu(cmn, 0) > 1.f) {
      kflfc += fem::isign(2, kflfc);
    }
  }
  else {
    kflfc = kflf(1);
  }
  lukfdi(cmn, kflfc, kflf(2), kfldmp, kf);
  if (kf == 0) {
    goto statement_280;
  }
  FEM_DO_SAFE(j, 1, fem::max(2, nfet)) {
    nfl(fem::iabs(kflf(j))) = nfl(fem::iabs(kflf(j))) - fem::isign(1, kflf(j));
  }
  //C
  //C...Store hadron at random among free positions.
  npos = fem::min(1 + fem::fint(rlu(cmn, 0) * nrem), nrem);
  FEM_DO_SAFE(i, nsav + njet + 1, n) {
    if (k(i, 1) == 7) {
      npos = npos - 1;
    }
    if (k(i, 1) == 1 || npos != 0) {
      goto statement_310;
    }
    k(i, 1) = 1;
    k(i, 2) = kf;
    p(i, 5) = ulmass(cmn, k(i, 2));
    p(i, 4) = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) +
      fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
    statement_310:;
  }
  nrem = nrem - 1;
  nreq = (fem::iabs(nfl(1)) + fem::iabs(nfl(2)) + fem::iabs(nfl(3)) -
    fem::iabs(nfl(1) + nfl(2) + nfl(3))) / 2 + fem::iabs(nfl(1) + nfl(2) +
    nfl(3)) / 3;
  if (nrem > 0) {
    goto statement_280;
  }
  //C
  //C...Compensate for missing momentum in global scheme (3 options).
  statement_320:
  if (fem::mod(mstj(3), 5) != 0 && fem::mod(mstj(3), 5) != 4) {
    FEM_DO_SAFE(j, 1, 3) {
      psi(j) = 0.f;
      FEM_DO_SAFE(i, nsav + njet + 1, n) {
        psi(j) += p(i, j);
      }
    }
    psi(4) = fem::pow2(psi(1)) + fem::pow2(psi(2)) + fem::pow2(psi(3));
    pws = 0.f;
    FEM_DO_SAFE(i, nsav + njet + 1, n) {
      if (fem::mod(mstj(3), 5) == 1) {
        pws += p(i, 4);
      }
      if (fem::mod(mstj(3), 5) == 2) {
        pws += fem::sqrt(fem::pow2(p(i, 5)) + fem::pow2((psi(1) * p(i,
          1) + psi(2) * p(i, 2) + psi(3) * p(i, 3))) / psi(4));
      }
      if (fem::mod(mstj(3), 5) == 3) {
        pws += 1.f;
      }
    }
    FEM_DO_SAFE(i, nsav + njet + 1, n) {
      if (fem::mod(mstj(3), 5) == 1) {
        pw = p(i, 4);
      }
      if (fem::mod(mstj(3), 5) == 2) {
        pw = fem::sqrt(fem::pow2(p(i, 5)) + fem::pow2((psi(1) * p(i,
          1) + psi(2) * p(i, 2) + psi(3) * p(i, 3))) / psi(4));
      }
      if (fem::mod(mstj(3), 5) == 3) {
        pw = 1.f;
      }
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) = p(i, j) - psi(j) * pw / pws;
      }
      p(i, 4) = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) +
        fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
    }
    //C
    //C...Compensate for missing momentum withing each jet separately.
  }
  else if (fem::mod(mstj(3), 5) == 4) {
    FEM_DO_SAFE(i, n + 1, n + njet) {
      k(i, 1) = 0;
      FEM_DO_SAFE(j, 1, 5) {
        p(i, j) = 0.f;
      }
    }
    FEM_DO_SAFE(i, nsav + njet + 1, n) {
      ir1 = k(i, 3);
      ir2 = n + ir1 - nsav;
      k(ir2, 1)++;
      pls = (p(i, 1) * p(ir1, 1) + p(i, 2) * p(ir1, 2) + p(i, 3) * p(ir1,
        3)) / (fem::pow2(p(ir1, 1)) + fem::pow2(p(ir1, 2)) + fem::pow2(p(ir1,
        3)));
      FEM_DO_SAFE(j, 1, 3) {
        p(ir2, j) += p(i, j) - pls * p(ir1, j);
      }
      p(ir2, 4) += p(i, 4);
      p(ir2, 5) += pls;
    }
    pss = 0.f;
    FEM_DO_SAFE(i, n + 1, n + njet) {
      if (k(i, 1) != 0) {
        pss += p(i, 4) / (pecm * (0.8f * p(i, 5) + 0.2f));
      }
    }
    FEM_DO_SAFE(i, nsav + njet + 1, n) {
      ir1 = k(i, 3);
      ir2 = n + ir1 - nsav;
      pls = (p(i, 1) * p(ir1, 1) + p(i, 2) * p(ir1, 2) + p(i, 3) * p(ir1,
        3)) / (fem::pow2(p(ir1, 1)) + fem::pow2(p(ir1, 2)) + fem::pow2(p(ir1,
        3)));
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) = p(i, j) - p(ir2, j) / k(ir2, 1) + (1.f / (p(ir2,
          5) * pss) - 1.f) * pls * p(ir1, j);
      }
      p(i, 4) = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) +
        fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
    }
  }
  //C
  //C...Scale momenta for energy conservation.
  if (fem::mod(mstj(3), 5) != 0) {
    pms = 0.f;
    pes = 0.f;
    pqs = 0.f;
    FEM_DO_SAFE(i, nsav + njet + 1, n) {
      pms += p(i, 5);
      pes += p(i, 4);
      pqs += fem::pow2(p(i, 5)) / p(i, 4);
    }
    if (pms >= pecm) {
      goto statement_150;
    }
    neco = 0;
    statement_440:
    neco++;
    pfac = (pecm - pqs) / (pes - pqs);
    pes = 0.f;
    pqs = 0.f;
    FEM_DO_SAFE(i, nsav + njet + 1, n) {
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) = pfac * p(i, j);
      }
      p(i, 4) = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) +
        fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
      pes += p(i, 4);
      pqs += fem::pow2(p(i, 5)) / p(i, 4);
    }
    if (neco < 10 && fem::abs(pecm - pes) > 2e-6f * pecm) {
      goto statement_440;
    }
  }
  //C
  //C...Origin of produced particles and parton daughter pointers.
  statement_470:
  FEM_DO_SAFE(i, nsav + njet + 1, n) {
    if (mstu(16) != 2) {
      k(i, 3) = nsav + 1;
    }
    if (mstu(16) == 2) {
      k(i, 3) = k(k(i, 3), 3);
    }
  }
  FEM_DO_SAFE(i, nsav + 1, nsav + njet) {
    i1 = k(i, 3);
    k(i1, 1) += 10;
    if (mstu(16) != 2) {
      k(i1, 4) = nsav + 1;
      k(i1, 5) = nsav + 1;
    }
    else {
      k(i1, 4) = k(i1, 4) - njet + 1;
      k(i1, 5) = k(i1, 5) - njet + 1;
      if (k(i1, 5) < k(i1, 4)) {
        k(i1, 4) = 0;
        k(i1, 5) = 0;
      }
    }
  }
  //C
  //C...Document independent fragmentation system. Remove copy of jets.
  nsav++;
  k(nsav, 1) = 11;
  k(nsav, 2) = 93;
  k(nsav, 3) = ip;
  k(nsav, 4) = nsav + 1;
  k(nsav, 5) = n - njet + 1;
  FEM_DO_SAFE(j, 1, 4) {
    p(nsav, j) = fem::sngl(dps(j));
    v(nsav, j) = v(ip, j);
  }
  p(nsav, 5) = fem::sqrt(fem::sngl(fem::max(0e0, fem::pow2(dps(4)) -
    fem::pow2(dps(1)) - fem::pow2(dps(2)) - fem::pow2(dps(3)))));
  v(nsav, 5) = 0.f;
  FEM_DO_SAFE(i, nsav + njet, n) {
    FEM_DO_SAFE(j, 1, 5) {
      k(i - njet + 1, j) = k(i, j);
      p(i - njet + 1, j) = p(i, j);
      v(i - njet + 1, j) = v(i, j);
    }
  }
  n = n - njet + 1;
  //C
  //C...Boost back particle system. Set production vertices.
  if (njet != 1) {
    ludbrb(nsav + 1, n, 0.f, 0.f, dps(1) / dps(4), dps(2) / dps(4),
      dps(3) / dps(4));
  }
  FEM_DO_SAFE(i, nsav + 1, n) {
    FEM_DO_SAFE(j, 1, 4) {
      v(i, j) = v(ip, j);
    }
  }
  //C
}

struct ludecy_save
{
  arr<float> wtcor;

  ludecy_save() :
    wtcor(dimension(10), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
ludecy(
  common& cmn,
  int const& ip)
{
  FEM_CMN_SVE(ludecy);
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_cref<float> brat(cmn.brat, dimension(2000));
  arr_cref<int, 2> kfdp(cmn.kfdp, dimension(2000, 5));
  int& nsav = cmn.nsav;
  //
  arr_ref<float> wtcor(sve.wtcor, dimension(10));
  if (is_called_first_time) {
    static const float values[] = {
      2.f, 5.f, 15.f, 60.f, 250.f, 1500.f, 1.2e4f, 1.2e5f, 150.f, 16.f
    };
    fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
      wtcor;
  }
  float a = fem::float0;
  float b = fem::float0;
  float c = fem::float0;
  int i = fem::int0;
  int j = fem::int0;
  float ha = fem::float0;
  float hrq = fem::float0;
  int ntry = fem::int0;
  int kfa = fem::int0;
  int kfs = fem::int0;
  int kc = fem::int0;
  arr_1d<4, float> vdcy(fem::fill0);
  int mout = fem::int0;
  int kca = fem::int0;
  int mdmdcy = fem::int0;
  int kfsp = fem::int0;
  int kfsn = fem::int0;
  int nope = fem::int0;
  float brsu = fem::float0;
  int idl = fem::int0;
  float rbr = fem::float0;
  int idc = fem::int0;
  int mmat = fem::int0;
  int np = fem::int0;
  int nq = fem::int0;
  int mbst = fem::int0;
  arr_2d<10, 5, float> pv(fem::fill0);
  float ps = fem::float0;
  float psq = fem::float0;
  int mrem = fem::int0;
  int jtmax = fem::int0;
  int jt = fem::int0;
  int kp = fem::int0;
  int kpa = fem::int0;
  int kcp = fem::int0;
  int kfp = fem::int0;
  int kdump = fem::int0;
  int kfpa = fem::int0;
  int kqp = fem::int0;
  arr_1d<4, int> kflo(fem::fill0);
  int kfi = fem::int0;
  int kfldmp = fem::int0;
  float psp = fem::float0;
  float cnde = fem::float0;
  float gauss = fem::float0;
  int nd = fem::int0;
  arr_1d<4, int> kfl1(fem::fill0);
  int kfl2 = fem::int0;
  int jt2 = fem::int0;
  int jt3 = fem::int0;
  float pqt = fem::float0;
  float psmc = fem::float0;
  float hr1 = fem::float0;
  float hr2 = fem::float0;
  float hlq = fem::float0;
  float huq = fem::float0;
  float hw = fem::float0;
  float hqw = fem::float0;
  float hlw = fem::float0;
  float huw = fem::float0;
  float hg = fem::float0;
  float hatl = fem::float0;
  float hm = fem::float0;
  float hmv1 = fem::float0;
  float hmv2 = fem::float0;
  float hsav1 = fem::float0;
  float hsav2 = fem::float0;
  float hmv = fem::float0;
  float hm1 = fem::float0;
  float hatm = fem::float0;
  float hwt1 = fem::float0;
  float hwt2 = fem::float0;
  float hwt3 = fem::float0;
  float hatu = fem::float0;
  float hmp1 = fem::float0;
  float hreg = fem::float0;
  float hacc = fem::float0;
  int nm = fem::int0;
  int msgn = fem::int0;
  int im = fem::int0;
  int kfam = fem::int0;
  int il = fem::int0;
  float wtmax = fem::float0;
  float pmax = fem::float0;
  float pmin = fem::float0;
  float pmes = fem::float0;
  float pmrho2 = fem::float0;
  float pgrho2 = fem::float0;
  float pmst = fem::float0;
  float wt = fem::float0;
  arr_1d<10, float> rord(fem::fill0);
  int il1 = fem::int0;
  float rsav = fem::float0;
  int il2 = fem::int0;
  float pa = fem::float0;
  arr_1d<3, float> ue(fem::fill0);
  float phi = fem::float0;
  arr_1d<3, float> be(fem::fill0);
  float ga = fem::float0;
  float bep = fem::float0;
  float four12 = fem::float0;
  float four13 = fem::float0;
  float four23 = fem::float0;
  float hx1 = fem::float0;
  float hx2 = fem::float0;
  float hx3 = fem::float0;
  int is = fem::int0;
  float pm2 = fem::float0;
  float pm3 = fem::float0;
  int kftemp = fem::int0;
  float pm4 = fem::float0;
  float hb = fem::float0;
  float hc = fem::float0;
  float hd = fem::float0;
  float he = fem::float0;
  float hf = fem::float0;
  float hh = fem::float0;
  float pcor = fem::float0;
  float pmr = fem::float0;
  float pm1 = fem::float0;
  int kfldum = fem::int0;
  int kf1 = fem::int0;
  int kf2 = fem::int0;
  float psm = fem::float0;
  int jcon = fem::int0;
  //C
  //C...Purpose: to handle the decay of unstable particles.
  //Clin-2/18/03 for resonance decay in hadron cascade:
  //C
  //C...Functions: momentum in two-particle decays, four-product and
  //C...matrix element times phase space in weak decays.
  pawt(a, b, c) = fem::sqrt((fem::pow2(a) - fem::pow2((b + c))) * (
    fem::pow2(a) - fem::pow2((b - c)))) / (2.f * a);
  four(i, j) = p(i, 4) * p(j, 4) - p(i, 1) * p(j, 1) - p(i, 2) * p(j,
    2) - p(i, 3) * p(j, 3);
  hmeps(ha) = (fem::pow2((1.f - hrq - ha)) + 3.f * ha * (1.f + hrq -
    ha)) * fem::sqrt(fem::pow2((1.f - hrq - ha)) - 4.f * hrq * ha);
  //C
  //C...Initial values.
  ntry = 0;
  nsav = n;
  kfa = fem::iabs(k(ip, 2));
  kfs = fem::isign(1, k(ip, 2));
  kc = lucomp(cmn, kfa);
  mstj(92) = 0;
  //C
  //C...Choose lifetime and determine decay vertex.
  if (k(ip, 1) == 5) {
    v(ip, 5) = 0.f;
  }
  else if (k(ip, 1) != 4) {
    v(ip, 5) = -pmas(kc, 4) * fem::log(rlu(cmn, 0));
  }
  FEM_DO_SAFE(j, 1, 4) {
    vdcy(j) = v(ip, j) + v(ip, 5) * p(ip, j) / p(ip, 5);
  }
  //C
  //C...Determine whether decay allowed or not.
  mout = 0;
  if (mstj(22) == 2) {
    if (pmas(kc, 4) > parj(71)) {
      mout = 1;
    }
  }
  else if (mstj(22) == 3) {
    if (fem::pow2(vdcy(1)) + fem::pow2(vdcy(2)) + fem::pow2(vdcy(
        3)) > fem::pow2(parj(72))) {
      mout = 1;
    }
  }
  else if (mstj(22) == 4) {
    if (fem::pow2(vdcy(1)) + fem::pow2(vdcy(2)) > fem::pow2(parj(73))) {
      mout = 1;
    }
    if (fem::abs(vdcy(3)) > parj(74)) {
      mout = 1;
    }
  }
  if (mout == 1 && k(ip, 1) != 5) {
    k(ip, 1) = 4;
    return;
  }
  //C
  //C...Check existence of decay channels. Particle/antiparticle rules.
  kca = kc;
  if (mdcy(kc, 2) > 0) {
    mdmdcy = mdme(mdcy(kc, 2), 2);
    if (mdmdcy > 80 && mdmdcy <= 90) {
      kca = mdmdcy;
    }
  }
  if (mdcy(kca, 2) <= 0 || mdcy(kca, 3) <= 0) {
    luerrm(cmn, 9, "(LUDECY:) no decay channel defined");
    return;
  }
  if (fem::mod(kfa / 1000, 10) == 0 && (kca == 85 || kca == 87)) {
    kfs = -kfs;
  }
  if (kchg(kc, 3) == 0) {
    kfsp = 1;
    kfsn = 0;
    if (rlu(cmn, 0) > 0.5f) {
      kfs = -kfs;
    }
  }
  else if (kfs > 0) {
    kfsp = 1;
    kfsn = 0;
  }
  else {
    kfsp = 0;
    kfsn = 1;
  }
  //C
  //C...Sum branching ratios of allowed decay channels.
  //Clin  110 NOPE=0
  nope = 0;
  brsu = 0.f;
  FEM_DO_SAFE(idl, mdcy(kca, 2), mdcy(kca, 2) + mdcy(kca, 3) - 1) {
    if (mdme(idl, 1) != 1 && kfsp * mdme(idl, 1) != 2 && kfsn * mdme(idl,
        1) != 3) {
      goto statement_120;
    }
    if (mdme(idl, 2) > 100) {
      goto statement_120;
    }
    nope++;
    brsu += brat(idl);
    statement_120:;
  }
  if (nope == 0) {
    luerrm(cmn, 2, "(LUDECY:) all decay channels closed by user");
    return;
  }
  //C
  //C...Select decay channel among allowed ones.
  statement_130:
  rbr = brsu * rlu(cmn, 0);
  idl = mdcy(kca, 2) - 1;
  statement_140:
  idl++;
  if (mdme(idl, 1) != 1 && kfsp * mdme(idl, 1) != 2 && kfsn * mdme(idl,
      1) != 3) {
    if (idl < mdcy(kca, 2) + mdcy(kca, 3) - 1) {
      goto statement_140;
    }
  }
  else if (mdme(idl, 2) > 100) {
    if (idl < mdcy(kca, 2) + mdcy(kca, 3) - 1) {
      goto statement_140;
    }
  }
  else {
    idc = idl;
    rbr = rbr - brat(idl);
    if (idl < mdcy(kca, 2) + mdcy(kca, 3) - 1 && rbr > 0.f) {
      goto statement_140;
    }
  }
  //C
  //C...Start readout of decay channel: matrix element, reset counters.
  mmat = mdme(idc, 2);
  statement_150:
  ntry++;
  if (ntry > 1000) {
    luerrm(cmn, 14, "(LUDECY:) caught in infinite loop");
    if (mstu(21) >= 1) {
      return;
    }
  }
  i = n;
  np = 0;
  nq = 0;
  mbst = 0;
  if (mmat >= 11 && mmat != 46 && p(ip, 4) > 20.f * p(ip, 5)) {
    mbst = 1;
  }
  FEM_DO_SAFE(j, 1, 4) {
    pv(1, j) = 0.f;
    if (mbst == 0) {
      pv(1, j) = p(ip, j);
    }
  }
  if (mbst == 1) {
    pv(1, 4) = p(ip, 5);
  }
  pv(1, 5) = p(ip, 5);
  ps = 0.f;
  psq = 0.f;
  mrem = 0;
  //C
  //C...Read out decay products. Convert to standard flavour code.
  jtmax = 5;
  if (mdme(idc + 1, 2) == 101) {
    jtmax = 10;
  }
  FEM_DO_SAFE(jt, 1, jtmax) {
    if (jt <= 5) {
      kp = kfdp(idc, jt);
    }
    if (jt >= 6) {
      kp = kfdp(idc + 1, jt - 5);
    }
    if (kp == 0) {
      goto statement_170;
    }
    kpa = fem::iabs(kp);
    kcp = lucomp(cmn, kpa);
    if (kchg(kcp, 3) == 0 && kpa != 81 && kpa != 82) {
      kfp = kp;
    }
    else if (kpa != 81 && kpa != 82) {
      kfp = kfs * kp;
    }
    else if (kpa == 81 && fem::mod(kfa / 1000, 10) == 0) {
      kfp = -kfs * fem::mod(kfa / 10, 10);
    }
    else if (kpa == 81 && fem::mod(kfa / 100, 10) >= fem::mod(kfa / 10, 10)) {
      kfp = kfs * (100 * fem::mod(kfa / 10, 100) + 3);
    }
    else if (kpa == 81) {
      kfp = kfs * (1000 * fem::mod(kfa / 10, 10) + 100 * fem::mod(kfa / 100,
        10) + 1);
    }
    else if (kp == 82) {
      lukfdi(cmn, -kfs * fem::fint(1.f + (2.f + parj(2)) * rlu(cmn,
        0)), 0, kfp, kdump);
      if (kfp == 0) {
        goto statement_150;
      }
      mstj(93) = 1;
      if (pv(1, 5) < parj(32) + 2.f * ulmass(cmn, kfp)) {
        goto statement_150;
      }
    }
    else if (kp ==  - 82) {
      kfp = -kfp;
      if (fem::iabs(kfp) > 10) {
        kfp += fem::isign(10000, kfp);
      }
    }
    if (kpa == 81 || kpa == 82) {
      kcp = lucomp(cmn, kfp);
    }
    //C
    //C...Add decay product to event record or to quark flavour list.
    kfpa = fem::iabs(kfp);
    kqp = kchg(kcp, 2);
    if (mmat >= 11 && mmat <= 30 && kqp != 0) {
      nq++;
      kflo(nq) = kfp;
      mstj(93) = 2;
      psq += ulmass(cmn, kflo(nq));
    }
    else if (mmat >= 42 && mmat <= 43 && np == 3 && fem::mod(nq, 2) == 1) {
      nq = nq - 1;
      ps = ps - p(i, 5);
      k(i, 1) = 1;
      kfi = k(i, 2);
      lukfdi(cmn, kfp, kfi, kfldmp, k(i, 2));
      if (k(i, 2) == 0) {
        goto statement_150;
      }
      mstj(93) = 1;
      p(i, 5) = ulmass(cmn, k(i, 2));
      ps += p(i, 5);
    }
    else {
      i++;
      np++;
      if (mmat != 33 && kqp != 0) {
        nq++;
      }
      if (mmat == 33 && kqp != 0 && kqp != 2) {
        nq++;
      }
      k(i, 1) = 1 + fem::mod(nq, 2);
      if (mmat == 4 && jt <= 2 && kfp == 21) {
        k(i, 1) = 2;
      }
      if (mmat == 4 && jt == 3) {
        k(i, 1) = 1;
      }
      k(i, 2) = kfp;
      k(i, 3) = ip;
      k(i, 4) = 0;
      k(i, 5) = 0;
      p(i, 5) = ulmass(cmn, kfp);
      if (mmat == 45 && kfpa == 89) {
        p(i, 5) = parj(32);
      }
      ps += p(i, 5);
    }
    statement_170:;
  }
  //C
  //C...Choose decay multiplicity in phase space model.
  statement_180:
  if (mmat >= 11 && mmat <= 30) {
    psp = ps;
    cnde = parj(61) * fem::log(fem::max((pv(1, 5) - ps - psq) / parj(62),
      1.1f));
    if (mmat == 12) {
      cnde += parj(63);
    }
    statement_190:
    ntry++;
    if (ntry > 1000) {
      luerrm(cmn, 14, "(LUDECY:) caught in infinite loop");
      if (mstu(21) >= 1) {
        return;
      }
    }
    if (mmat <= 20) {
      gauss = fem::sqrt(-2.f * cnde * fem::log(fem::max(1e-10f, rlu(cmn,
        0)))) * fem::sin(paru(2) * rlu(cmn, 0));
      nd = fem::fint(0.5f + 0.5f * np + 0.25f * nq + cnde + gauss);
      if (nd < np + nq / 2 || nd < 2 || nd > 10) {
        goto statement_190;
      }
      if (mmat == 13 && nd == 2) {
        goto statement_190;
      }
      if (mmat == 14 && nd <= 3) {
        goto statement_190;
      }
      if (mmat == 15 && nd <= 4) {
        goto statement_190;
      }
    }
    else {
      nd = mmat - 20;
    }
    //C
    //C...Form hadrons from flavour content.
    FEM_DO_SAFE(jt, 1, 4) {
      kfl1(jt) = kflo(jt);
    }
    if (nd == np + nq / 2) {
      goto statement_220;
    }
    FEM_DO_SAFE(i, n + np + 1, n + nd - nq / 2) {
      jt = 1 + fem::fint((nq - 1) * rlu(cmn, 0));
      lukfdi(cmn, kfl1(jt), 0, kfl2, k(i, 2));
      if (k(i, 2) == 0) {
        goto statement_190;
      }
      kfl1(jt) = -kfl2;
    }
    statement_220:
    jt = 2;
    jt2 = 3;
    jt3 = 4;
    if (nq == 4 && rlu(cmn, 0) < parj(66)) {
      jt = 4;
    }
    if (jt == 4 && fem::isign(1, kfl1(1) * (10 - fem::iabs(kfl1(
        1)))) * fem::isign(1, kfl1(jt) * (10 - fem::iabs(kfl1(
        jt)))) > 0) {
      jt = 3;
    }
    if (jt == 3) {
      jt2 = 2;
    }
    if (jt == 4) {
      jt3 = 2;
    }
    lukfdi(cmn, kfl1(1), kfl1(jt), kfldmp, k(n + nd - nq / 2 + 1, 2));
    if (k(n + nd - nq / 2 + 1, 2) == 0) {
      goto statement_190;
    }
    if (nq == 4) {
      lukfdi(cmn, kfl1(jt2), kfl1(jt3), kfldmp, k(n + nd, 2));
    }
    if (nq == 4 && k(n + nd, 2) == 0) {
      goto statement_190;
    }
    //C
    //C...Check that sum of decay product masses not too large.
    ps = psp;
    FEM_DO_SAFE(i, n + np + 1, n + nd) {
      k(i, 1) = 1;
      k(i, 3) = ip;
      k(i, 4) = 0;
      k(i, 5) = 0;
      p(i, 5) = ulmass(cmn, k(i, 2));
      ps += p(i, 5);
    }
    if (ps + parj(64) > pv(1, 5)) {
      goto statement_190;
    }
    //C
    //C...Rescale energy to subtract off spectator quark mass.
  }
  else if ((mmat == 31 || mmat == 33 || mmat == 44 || mmat == 45) && np >= 3) {
    ps = ps - p(n + np, 5);
    pqt = (p(n + np, 5) + parj(65)) / pv(1, 5);
    FEM_DO_SAFE(j, 1, 5) {
      p(n + np, j) = pqt * pv(1, j);
      pv(1, j) = (1.f - pqt) * pv(1, j);
    }
    if (ps + parj(64) > pv(1, 5)) {
      goto statement_150;
    }
    nd = np - 1;
    mrem = 1;
    //C
    //C...Phase space factors imposed in W decay.
  }
  else if (mmat == 46) {
    mstj(93) = 1;
    psmc = ulmass(cmn, k(n + 1, 2));
    mstj(93) = 1;
    psmc += ulmass(cmn, k(n + 2, 2));
    if (fem::max(ps, psmc) + parj(32) > pv(1, 5)) {
      goto statement_130;
    }
    hr1 = fem::pow2((p(n + 1, 5) / pv(1, 5)));
    hr2 = fem::pow2((p(n + 2, 5) / pv(1, 5)));
    if ((1.f - hr1 - hr2) * (2.f + hr1 + hr2) * fem::sqrt(fem::pow2((
        1.f - hr1 - hr2)) - 4.f * hr1 * hr2) < 2.f * rlu(cmn, 0)) {
      goto statement_130;
    }
    nd = np;
    //C
    //C...Fully specified final state: check mass broadening effects.
  }
  else {
    if (np >= 2 && ps + parj(64) > pv(1, 5)) {
      goto statement_150;
    }
    nd = np;
  }
  //C
  //C...Select W mass in decay Q -> W + q, without W propagator.
  if (mmat == 45 && mstj(25) <= 0) {
    hlq = fem::pow2((parj(32) / pv(1, 5)));
    huq = fem::pow2((1.f - (p(n + 2, 5) + parj(64)) / pv(1, 5)));
    hrq = fem::pow2((p(n + 2, 5) / pv(1, 5)));
    statement_250:
    hw = hlq + rlu(cmn, 0) * (huq - hlq);
    if (hmeps(hw) < rlu(cmn, 0)) {
      goto statement_250;
    }
    p(n + 1, 5) = pv(1, 5) * fem::sqrt(hw);
    //C
    //C...Ditto, including W propagator. Divide mass range into three regions.
  }
  else if (mmat == 45) {
    hqw = fem::pow2((pv(1, 5) / pmas(24, 1)));
    hlw = fem::pow2((parj(32) / pmas(24, 1)));
    huw = fem::pow2(((pv(1, 5) - p(n + 2, 5) - parj(64)) / pmas(24, 1)));
    hrq = fem::pow2((p(n + 2, 5) / pv(1, 5)));
    hg = pmas(24, 2) / pmas(24, 1);
    hatl = fem::atan((hlw - 1.f) / hg);
    hm = fem::min(1.f, huw - 0.001f);
    hmv1 = hmeps(hm / hqw) / (fem::pow2((hm - 1.f)) + fem::pow2(hg));
    statement_260:
    hm = hm - hg;
    hmv2 = hmeps(hm / hqw) / (fem::pow2((hm - 1.f)) + fem::pow2(hg));
    hsav1 = hmeps(hm / hqw);
    hsav2 = 1.f / (fem::pow2((hm - 1.f)) + fem::pow2(hg));
    if (hmv2 > hmv1 && hm - hg > hlw) {
      hmv1 = hmv2;
      goto statement_260;
    }
    hmv = fem::min(2.f * hmv1, hmeps(hm / hqw) / fem::pow2(hg));
    hm1 = 1.f - fem::sqrt(1.f / hmv - fem::pow2(hg));
    if (hm1 > hlw && hm1 < hm) {
      hm = hm1;
    }
    else if (hmv2 <= hmv1) {
      hm = fem::max(hlw, hm - fem::min(0.1f, 1.f - hm));
    }
    hatm = fem::atan((hm - 1.f) / hg);
    hwt1 = (hatm - hatl) / hg;
    hwt2 = hmv * (fem::min(1.f, huw) - hm);
    hwt3 = 0.f;
    if (huw > 1.f) {
      hatu = fem::atan((huw - 1.f) / hg);
      hmp1 = hmeps(1.f / hqw);
      hwt3 = hmp1 * hatu / hg;
    }
    //C
    //C...Select mass region and W mass there. Accept according to weight.
    statement_270:
    hreg = rlu(cmn, 0) * (hwt1 + hwt2 + hwt3);
    if (hreg <= hwt1) {
      hw = 1.f + hg * fem::tan(hatl + rlu(cmn, 0) * (hatm - hatl));
      hacc = hmeps(hw / hqw);
    }
    else if (hreg <= hwt1 + hwt2) {
      hw = hm + rlu(cmn, 0) * (fem::min(1.f, huw) - hm);
      hacc = hmeps(hw / hqw) / (fem::pow2((hw - 1.f)) + fem::pow2(hg)) / hmv;
    }
    else {
      hw = 1.f + hg * fem::tan(rlu(cmn, 0) * hatu);
      hacc = hmeps(hw / hqw) / hmp1;
    }
    if (hacc < rlu(cmn, 0)) {
      goto statement_270;
    }
    p(n + 1, 5) = pmas(24, 1) * fem::sqrt(hw);
  }
  //C
  //C...Determine position of grandmother, number of sisters, Q -> W sign.
  nm = 0;
  msgn = 0;
  if (mmat == 3 || mmat == 46) {
    im = k(ip, 3);
    if (im < 0 || im >= ip) {
      im = 0;
    }
    if (im != 0) {
      kfam = fem::iabs(k(im, 2));
    }
    if (im != 0 && mmat == 3) {
      FEM_DO_SAFE(il, fem::max(ip - 2, im + 1), fem::min(ip + 2, n)) {
        if (k(il, 3) == im) {
          nm++;
        }
      }
      if (nm != 2 || kfam <= 100 || fem::mod(kfam, 10) != 1 ||
          fem::mod(kfam / 1000, 10) != 0) {
        nm = 0;
      }
    }
    else if (im != 0 && mmat == 46) {
      msgn = fem::isign(1, k(im, 2) * k(ip, 2));
      if (kfam > 100 && fem::mod(kfam / 1000, 10) == 0) {
        msgn = msgn * fem::pow((-1), fem::mod(kfam / 100, 10));
      }
    }
  }
  //C
  //C...Kinematics of one-particle decays.
  if (nd == 1) {
    FEM_DO_SAFE(j, 1, 4) {
      p(n + 1, j) = p(ip, j);
    }
    goto statement_510;
  }
  //C
  //C...Calculate maximum weight ND-particle decay.
  pv(nd, 5) = p(n + nd, 5);
  if (nd >= 3) {
    wtmax = 1.f / wtcor(nd - 2);
    pmax = pv(1, 5) - ps + p(n + nd, 5);
    pmin = 0.f;
    FEM_DOSTEP(il, nd - 1, 1, -1) {
      pmax += p(n + il, 5);
      pmin += p(n + il + 1, 5);
      wtmax = wtmax * pawt(pmax, pmin, p(n + il, 5));
    }
  }
  //C
  //C...Find virtual gamma mass in Dalitz decay.
  statement_310:
  if (nd == 2) {
  }
  else if (mmat == 2) {
    pmes = 4.f * fem::pow2(pmas(11, 1));
    pmrho2 = fem::pow2(pmas(131, 1));
    pgrho2 = fem::pow2(pmas(131, 2));
    statement_320:
    pmst = pmes * fem::pow((fem::pow2(p(ip, 5)) / pmes), rlu(cmn, 0));
    wt = (1 + 0.5f * pmes / pmst) * fem::sqrt(fem::max(0.f, 1.f -
      pmes / pmst)) * fem::pow3((1.f - pmst / fem::pow2(p(ip, 5)))) *
      (1.f + pgrho2 / pmrho2) / (fem::pow2((1.f - pmst / pmrho2)) +
      pgrho2 / pmrho2);
    if (wt < rlu(cmn, 0)) {
      goto statement_320;
    }
    pv(2, 5) = fem::max(2.00001f * pmas(11, 1), fem::sqrt(pmst));
    //C
    //C...M-generator gives weight. If rejected, try again.
  }
  else {
    statement_330:
    rord(1) = 1.f;
    FEM_DO_SAFE(il1, 2, nd - 1) {
      rsav = rlu(cmn, 0);
      FEM_DOSTEP(il2, il1 - 1, 1, -1) {
        if (rsav <= rord(il2)) {
          goto statement_350;
        }
        rord(il2 + 1) = rord(il2);
      }
      statement_350:
      rord(il2 + 1) = rsav;
    }
    rord(nd) = 0.f;
    wt = 1.f;
    FEM_DOSTEP(il, nd - 1, 1, -1) {
      pv(il, 5) = pv(il + 1, 5) + p(n + il, 5) + (rord(il) - rord(
        il + 1)) * (pv(1, 5) - ps);
      wt = wt * pawt(pv(il, 5), pv(il + 1, 5), p(n + il, 5));
    }
    if (wt < rlu(cmn, 0) * wtmax) {
      goto statement_330;
    }
  }
  //C
  //C...Perform two-particle decays in respective CM frame.
  statement_370:
  FEM_DO_SAFE(il, 1, nd - 1) {
    pa = pawt(pv(il, 5), pv(il + 1, 5), p(n + il, 5));
    ue(3) = 2.f * rlu(cmn, 0) - 1.f;
    phi = paru(2) * rlu(cmn, 0);
    ue(1) = fem::sqrt(1.f - fem::pow2(ue(3))) * fem::cos(phi);
    ue(2) = fem::sqrt(1.f - fem::pow2(ue(3))) * fem::sin(phi);
    FEM_DO_SAFE(j, 1, 3) {
      p(n + il, j) = pa * ue(j);
      pv(il + 1, j) = -pa * ue(j);
    }
    p(n + il, 4) = fem::sqrt(fem::pow2(pa) + fem::pow2(p(n + il, 5)));
    pv(il + 1, 4) = fem::sqrt(fem::pow2(pa) + fem::pow2(pv(il + 1, 5)));
  }
  //C
  //C...Lorentz transform decay products to lab frame.
  FEM_DO_SAFE(j, 1, 4) {
    p(n + nd, j) = pv(nd, j);
  }
  FEM_DOSTEP(il, nd - 1, 1, -1) {
    FEM_DO_SAFE(j, 1, 3) {
      be(j) = pv(il, j) / pv(il, 4);
    }
    ga = pv(il, 4) / pv(il, 5);
    FEM_DO_SAFE(i, n + il, n + nd) {
      bep = be(1) * p(i, 1) + be(2) * p(i, 2) + be(3) * p(i, 3);
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) += ga * (ga * bep / (1.f + ga) + p(i, 4)) * be(j);
      }
      p(i, 4) = ga * (p(i, 4) + bep);
    }
  }
  //C
  //C...Matrix elements for omega and phi decays.
  if (mmat == 1) {
    wt = fem::pow2((p(n + 1, 5) * p(n + 2, 5) * p(n + 3, 5))) -
      fem::pow2((p(n + 1, 5) * four(n + 2, n + 3))) - fem::pow2((p(n + 2,
      5) * four(n + 1, n + 3))) - fem::pow2((p(n + 3, 5) * four(n + 1,
      n + 2))) + 2.f * four(n + 1, n + 2) * four(n + 1, n + 3) * four(n + 2,
      n + 3);
    if (fem::max(wt * wtcor(9) / fem::pow(p(ip, 5), 6), 0.001f) < rlu(cmn, 0)) {
      goto statement_310;
    }
    //C
    //C...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-.
  }
  else if (mmat == 2) {
    four12 = four(n + 1, n + 2);
    four13 = four(n + 1, n + 3);
    four23 = 0.5f * pmst - 0.25f * pmes;
    wt = (pmst - 0.5f * pmes) * (fem::pow2(four12) + fem::pow2(
      four13)) + pmes * (four12 * four13 + fem::pow2(four12) +
      fem::pow2(four13));
    if (wt < rlu(cmn, 0) * 0.25f * pmst * fem::pow2((fem::pow2(p(ip,
        5)) - pmst))) {
      goto statement_370;
    }
    //C
    //C...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar,
    //C...V vector), of form cos**2(theta02) in V1 rest frame.
  }
  else if (mmat == 3 && nm == 2) {
    if (fem::pow2((fem::pow2(p(ip, 5)) * four(im, n + 1) - four(ip,
        im) * four(ip, n + 1))) <= rlu(cmn, 0) * (fem::pow2(four(ip,
        im)) - fem::pow2((p(ip, 5) * p(im, 5)))) * (fem::pow2(four(ip,
        n + 1)) - fem::pow2((p(ip, 5) * p(n + 1, 5))))) {
      goto statement_370;
    }
    //C
    //C...Matrix element for "onium" -> g + g + g or gamma + g + g.
  }
  else if (mmat == 4) {
    hx1 = 2.f * four(ip, n + 1) / fem::pow2(p(ip, 5));
    hx2 = 2.f * four(ip, n + 2) / fem::pow2(p(ip, 5));
    hx3 = 2.f * four(ip, n + 3) / fem::pow2(p(ip, 5));
    wt = fem::pow2(((1.f - hx1) / (hx2 * hx3))) + fem::pow2(((1.f -
      hx2) / (hx1 * hx3))) + fem::pow2(((1.f - hx3) / (hx1 * hx2)));
    if (wt < 2.f * rlu(cmn, 0)) {
      goto statement_310;
    }
    if (k(ip + 1, 2) == 22 && (1.f - hx1) * fem::pow2(p(ip,
        5)) < 4.f * fem::pow2(parj(32))) {
      goto statement_310;
    }
    //C
    //C...Effective matrix element for nu spectrum in tau -> nu + hadrons.
  }
  else if (mmat == 41) {
    hx1 = 2.f * four(ip, n + 1) / fem::pow2(p(ip, 5));
    if (8.f * hx1 * (3.f - 2.f * hx1) / 9.f < rlu(cmn, 0)) {
      goto statement_310;
    }
    //C
    //C...Matrix elements for weak decays (only semileptonic for c and b)
  }
  else if (mmat >= 42 && mmat <= 44 && nd == 3) {
    if (mbst == 0) {
      wt = four(ip, n + 1) * four(n + 2, n + 3);
    }
    if (mbst == 1) {
      wt = p(ip, 5) * p(n + 1, 4) * four(n + 2, n + 3);
    }
    if (wt < rlu(cmn, 0) * p(ip, 5) * fem::pow3(pv(1, 5)) / wtcor(10)) {
      goto statement_310;
    }
  }
  else if (mmat >= 42 && mmat <= 44) {
    FEM_DO_SAFE(j, 1, 4) {
      p(n + np + 1, j) = 0.f;
      FEM_DO_SAFE(is, n + 3, n + np) {
        p(n + np + 1, j) += p(is, j);
      }
    }
    if (mbst == 0) {
      wt = four(ip, n + 1) * four(n + 2, n + np + 1);
    }
    if (mbst == 1) {
      wt = p(ip, 5) * p(n + 1, 4) * four(n + 2, n + np + 1);
    }
    if (wt < rlu(cmn, 0) * p(ip, 5) * fem::pow3(pv(1, 5)) / wtcor(10)) {
      goto statement_310;
    }
    //C
    //C...Angular distribution in W decay.
  }
  else if (mmat == 46 && msgn != 0) {
    if (msgn > 0) {
      wt = four(im, n + 1) * four(n + 2, ip + 1);
    }
    if (msgn < 0) {
      wt = four(im, n + 2) * four(n + 1, ip + 1);
    }
    if (wt < rlu(cmn, 0) * fem::pow4(p(im, 5)) / wtcor(10)) {
      goto statement_370;
    }
  }
  //C
  //C...Scale back energy and reattach spectator.
  if (mrem == 1) {
    FEM_DO_SAFE(j, 1, 5) {
      pv(1, j) = pv(1, j) / (1.f - pqt);
    }
    nd++;
    mrem = 0;
  }
  //C
  //C...Low invariant mass for system with spectator quark gives particle,
  //C...not two jets. Readjust momenta accordingly.
  if ((mmat == 31 || mmat == 45) && nd == 3) {
    mstj(93) = 1;
    pm2 = ulmass(cmn, k(n + 2, 2));
    mstj(93) = 1;
    pm3 = ulmass(cmn, k(n + 3, 2));
    if (fem::pow2(p(n + 2, 5)) + fem::pow2(p(n + 3, 5)) + 2.f * four(n + 2,
        n + 3) >= fem::pow2((parj(32) + pm2 + pm3))) {
      goto statement_510;
    }
    k(n + 2, 1) = 1;
    kftemp = k(n + 2, 2);
    lukfdi(cmn, kftemp, k(n + 3, 2), kfldmp, k(n + 2, 2));
    if (k(n + 2, 2) == 0) {
      goto statement_150;
    }
    p(n + 2, 5) = ulmass(cmn, k(n + 2, 2));
    ps = p(n + 1, 5) + p(n + 2, 5);
    pv(2, 5) = p(n + 2, 5);
    mmat = 0;
    nd = 2;
    goto statement_370;
  }
  else if (mmat == 44) {
    mstj(93) = 1;
    pm3 = ulmass(cmn, k(n + 3, 2));
    mstj(93) = 1;
    pm4 = ulmass(cmn, k(n + 4, 2));
    if (fem::pow2(p(n + 3, 5)) + fem::pow2(p(n + 4, 5)) + 2.f * four(n + 3,
        n + 4) >= fem::pow2((parj(32) + pm3 + pm4))) {
      goto statement_480;
    }
    k(n + 3, 1) = 1;
    kftemp = k(n + 3, 2);
    lukfdi(cmn, kftemp, k(n + 4, 2), kfldmp, k(n + 3, 2));
    if (k(n + 3, 2) == 0) {
      goto statement_150;
    }
    p(n + 3, 5) = ulmass(cmn, k(n + 3, 2));
    FEM_DO_SAFE(j, 1, 3) {
      p(n + 3, j) += p(n + 4, j);
    }
    p(n + 3, 4) = fem::sqrt(fem::pow2(p(n + 3, 1)) + fem::pow2(p(n + 3,
      2)) + fem::pow2(p(n + 3, 3)) + fem::pow2(p(n + 3, 5)));
    ha = fem::pow2(p(n + 1, 4)) - fem::pow2(p(n + 2, 4));
    hb = ha - (fem::pow2(p(n + 1, 5)) - fem::pow2(p(n + 2, 5)));
    hc = fem::pow2((p(n + 1, 1) - p(n + 2, 1))) + fem::pow2((p(n + 1,
      2) - p(n + 2, 2))) + fem::pow2((p(n + 1, 3) - p(n + 2, 3)));
    hd = fem::pow2((pv(1, 4) - p(n + 3, 4)));
    he = fem::pow2(ha) - 2.f * hd * (fem::pow2(p(n + 1, 4)) +
      fem::pow2(p(n + 2, 4))) + fem::pow2(hd);
    hf = hd * hc - fem::pow2(hb);
    hg = hd * hc - ha * hb;
    hh = (fem::sqrt(fem::pow2(hg) + he * hf) - hg) / (2.f * hf);
    FEM_DO_SAFE(j, 1, 3) {
      pcor = hh * (p(n + 1, j) - p(n + 2, j));
      p(n + 1, j) += pcor;
      p(n + 2, j) = p(n + 2, j) - pcor;
    }
    p(n + 1, 4) = fem::sqrt(fem::pow2(p(n + 1, 1)) + fem::pow2(p(n + 1,
      2)) + fem::pow2(p(n + 1, 3)) + fem::pow2(p(n + 1, 5)));
    p(n + 2, 4) = fem::sqrt(fem::pow2(p(n + 2, 1)) + fem::pow2(p(n + 2,
      2)) + fem::pow2(p(n + 2, 3)) + fem::pow2(p(n + 2, 5)));
    nd = nd - 1;
  }
  //C
  //C...Check invariant mass of W jets. May give one particle or start over.
  statement_480:
  if (mmat >= 42 && mmat <= 44 && fem::iabs(k(n + 1, 2)) < 10) {
    pmr = fem::sqrt(fem::max(0.f, fem::pow2(p(n + 1, 5)) + fem::pow2(p(n + 2,
      5)) + 2.f * four(n + 1, n + 2)));
    mstj(93) = 1;
    pm1 = ulmass(cmn, k(n + 1, 2));
    mstj(93) = 1;
    pm2 = ulmass(cmn, k(n + 2, 2));
    if (pmr > parj(32) + pm1 + pm2) {
      goto statement_490;
    }
    kfldum = fem::fint(1.5f + rlu(cmn, 0));
    lukfdi(cmn, k(n + 1, 2), -fem::isign(kfldum, k(n + 1, 2)), kfldmp, kf1);
    lukfdi(cmn, k(n + 2, 2), -fem::isign(kfldum, k(n + 2, 2)), kfldmp, kf2);
    if (kf1 == 0 || kf2 == 0) {
      goto statement_150;
    }
    psm = ulmass(cmn, kf1) + ulmass(cmn, kf2);
    if (mmat == 42 && pmr > parj(64) + psm) {
      goto statement_490;
    }
    if (mmat >= 43 && pmr > 0.2f * parj(32) + psm) {
      goto statement_490;
    }
    if (nd == 4 || kfa == 15) {
      goto statement_150;
    }
    k(n + 1, 1) = 1;
    kftemp = k(n + 1, 2);
    lukfdi(cmn, kftemp, k(n + 2, 2), kfldmp, k(n + 1, 2));
    if (k(n + 1, 2) == 0) {
      goto statement_150;
    }
    p(n + 1, 5) = ulmass(cmn, k(n + 1, 2));
    k(n + 2, 2) = k(n + 3, 2);
    p(n + 2, 5) = p(n + 3, 5);
    ps = p(n + 1, 5) + p(n + 2, 5);
    pv(2, 5) = p(n + 3, 5);
    mmat = 0;
    nd = 2;
    goto statement_370;
  }
  //C
  //C...Phase space decay of partons from W decay.
  statement_490:
  if (mmat == 42 && fem::iabs(k(n + 1, 2)) < 10) {
    kflo(1) = k(n + 1, 2);
    kflo(2) = k(n + 2, 2);
    k(n + 1, 1) = k(n + 3, 1);
    k(n + 1, 2) = k(n + 3, 2);
    FEM_DO_SAFE(j, 1, 5) {
      pv(1, j) = p(n + 1, j) + p(n + 2, j);
      p(n + 1, j) = p(n + 3, j);
    }
    pv(1, 5) = pmr;
    n++;
    np = 0;
    nq = 2;
    ps = 0.f;
    mstj(93) = 2;
    psq = ulmass(cmn, kflo(1));
    mstj(93) = 2;
    psq += ulmass(cmn, kflo(2));
    mmat = 11;
    goto statement_180;
  }
  //C
  //C...Boost back for rapidly moving particle.
  statement_510:
  n += nd;
  if (mbst == 1) {
    FEM_DO_SAFE(j, 1, 3) {
      be(j) = p(ip, j) / p(ip, 4);
    }
    ga = p(ip, 4) / p(ip, 5);
    FEM_DO_SAFE(i, nsav + 1, n) {
      bep = be(1) * p(i, 1) + be(2) * p(i, 2) + be(3) * p(i, 3);
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) += ga * (ga * bep / (1.f + ga) + p(i, 4)) * be(j);
      }
      p(i, 4) = ga * (p(i, 4) + bep);
    }
  }
  //C
  //C...Fill in position of decay vertex.
  FEM_DO_SAFE(i, nsav + 1, n) {
    FEM_DO_SAFE(j, 1, 4) {
      v(i, j) = vdcy(j);
    }
    v(i, 5) = 0.f;
  }
  //C
  //C...Set up for parton shower evolution from jets.
  if (mstj(23) >= 1 && mmat == 4 && k(nsav + 1, 2) == 21) {
    k(nsav + 1, 1) = 3;
    k(nsav + 2, 1) = 3;
    k(nsav + 3, 1) = 3;
    k(nsav + 1, 4) = mstu(5) * (nsav + 2);
    k(nsav + 1, 5) = mstu(5) * (nsav + 3);
    k(nsav + 2, 4) = mstu(5) * (nsav + 3);
    k(nsav + 2, 5) = mstu(5) * (nsav + 1);
    k(nsav + 3, 4) = mstu(5) * (nsav + 1);
    k(nsav + 3, 5) = mstu(5) * (nsav + 2);
    mstj(92) = -(nsav + 1);
  }
  else if (mstj(23) >= 1 && mmat == 4) {
    k(nsav + 2, 1) = 3;
    k(nsav + 3, 1) = 3;
    k(nsav + 2, 4) = mstu(5) * (nsav + 3);
    k(nsav + 2, 5) = mstu(5) * (nsav + 3);
    k(nsav + 3, 4) = mstu(5) * (nsav + 2);
    k(nsav + 3, 5) = mstu(5) * (nsav + 2);
    mstj(92) = nsav + 2;
  }
  else if (mstj(23) >= 1 && (mmat == 32 || mmat == 44 ||
    mmat == 46) && fem::iabs(k(nsav + 1, 2)) <= 10 && fem::iabs(k(nsav + 2,
    2)) <= 10) {
    k(nsav + 1, 1) = 3;
    k(nsav + 2, 1) = 3;
    k(nsav + 1, 4) = mstu(5) * (nsav + 2);
    k(nsav + 1, 5) = mstu(5) * (nsav + 2);
    k(nsav + 2, 4) = mstu(5) * (nsav + 1);
    k(nsav + 2, 5) = mstu(5) * (nsav + 1);
    mstj(92) = nsav + 1;
  }
  else if (mstj(23) >= 1 && mmat == 33 && fem::iabs(k(nsav + 2, 2)) == 21) {
    k(nsav + 1, 1) = 3;
    k(nsav + 2, 1) = 3;
    k(nsav + 3, 1) = 3;
    kcp = lucomp(cmn, k(nsav + 1, 2));
    kqp = kchg(kcp, 2) * fem::isign(1, k(nsav + 1, 2));
    jcon = 4;
    if (kqp < 0) {
      jcon = 5;
    }
    k(nsav + 1, jcon) = mstu(5) * (nsav + 2);
    k(nsav + 2, 9 - jcon) = mstu(5) * (nsav + 1);
    k(nsav + 2, jcon) = mstu(5) * (nsav + 3);
    k(nsav + 3, 9 - jcon) = mstu(5) * (nsav + 2);
    mstj(92) = nsav + 1;
  }
  else if (mstj(23) >= 1 && mmat == 33) {
    k(nsav + 1, 1) = 3;
    k(nsav + 3, 1) = 3;
    k(nsav + 1, 4) = mstu(5) * (nsav + 3);
    k(nsav + 1, 5) = mstu(5) * (nsav + 3);
    k(nsav + 3, 4) = mstu(5) * (nsav + 1);
    k(nsav + 3, 5) = mstu(5) * (nsav + 1);
    mstj(92) = nsav + 1;
  }
  //C
  //C...Mark decayed particle.
  if (k(ip, 1) == 5) {
    k(ip, 1) = 15;
  }
  if (k(ip, 1) <= 10) {
    k(ip, 1) = 11;
  }
  k(ip, 4) = nsav + 1;
  k(ip, 5) = n;
  //C
}

//C
//C*********************************************************************
//C
void
lushow(
  common& cmn,
  int const& ip1,
  int const& ip2,
  float const& qmax)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  arr_2d<5, 40, float> pmth(fem::fill0);
  float pmqth1 = fem::float0;
  float pmqth2 = fem::float0;
  int identifier_if = fem::int0;
  float pt2min = fem::float0;
  float alams = fem::float0;
  float alfm = fem::float0;
  int m3jc = fem::int0;
  int npa = fem::int0;
  arr_1d<4, int> ipa(fem::fill0);
  int i = fem::int0;
  int irej = fem::int0;
  int j = fem::int0;
  arr_1d<5, float> ps(fem::fill0);
  float pm = fem::float0;
  arr_1d<4, int> kfla(fem::fill0);
  arr_1d<4, float> pma(fem::fill0);
  int ns = fem::int0;
  int nep = fem::int0;
  int im = fem::int0;
  int kflm = fem::int0;
  int igm = fem::int0;
  int iau = fem::int0;
  int ip = fem::int0;
  arr_1d<4, int> kfld(fem::fill0);
  arr_1d<4, int> itry(fem::fill0);
  arr_1d<4, int> isl(fem::fill0);
  arr_1d<4, int> isi(fem::fill0);
  int islm = fem::int0;
  float pem = fem::float0;
  arr_1d<4, float> pmsd(fem::fill0);
  int inum = fem::int0;
  float rmax = fem::float0;
  float rpm = fem::float0;
  arr_1d<4, int> iep(fem::fill0);
  arr_1d<4, int> kfl(fem::fill0);
  float z = fem::float0;
  float pmed = fem::float0;
  float zc = fem::float0;
  float zce = fem::float0;
  float fbr = fem::float0;
  float fbre = fem::float0;
  float pms = fem::float0;
  float pm2 = fem::float0;
  float b0 = fem::float0;
  float pmsqcd = fem::float0;
  int mce = fem::int0;
  float pmsqed = fem::float0;
  int kflb = fem::int0;
  float pmq = fem::float0;
  float pmq0 = fem::float0;
  int kflgd1 = fem::int0;
  int kflgd2 = fem::int0;
  float ped = fem::float0;
  float pmqth3 = fem::float0;
  float pmq1 = fem::float0;
  float pmq2 = fem::float0;
  float zd = fem::float0;
  float zh = fem::float0;
  float zl = fem::float0;
  float zu = fem::float0;
  float x1 = fem::float0;
  float x2 = fem::float0;
  float x3 = fem::float0;
  int ki1 = fem::int0;
  int ki2 = fem::int0;
  float qf1 = fem::float0;
  float qf2 = fem::float0;
  float wshow = fem::float0;
  float wme = fem::float0;
  int maom = fem::int0;
  float zm = fem::float0;
  float the2id = fem::float0;
  int iaom = fem::int0;
  float the2im = fem::float0;
  float pa1s = fem::float0;
  float pa2s = fem::float0;
  float pa3s = fem::float0;
  float pts = fem::float0;
  int i1 = fem::int0;
  int kflda = fem::int0;
  int i2 = fem::int0;
  float pml = fem::float0;
  float zdr1 = fem::float0;
  float zdr2 = fem::float0;
  int mazip = fem::int0;
  int mazic = fem::int0;
  float ped1 = fem::float0;
  float pzm = fem::float0;
  float pmls = fem::float0;
  float pt = fem::float0;
  float hazip = fem::float0;
  float zau = fem::float0;
  float hazic = fem::float0;
  float zs = fem::float0;
  float zgm = fem::float0;
  float phi = fem::float0;
  float bex = fem::float0;
  float bey = fem::float0;
  float bez = fem::float0;
  float ga = fem::float0;
  float gabep = fem::float0;
  float the = fem::float0;
  arr_1d<4, double> dp(fem::fill0);
  double dbp = fem::double0;
  double dgabp = fem::double0;
  arr_2d<5, 4, double> dpt(fem::fill0);
  double dpma = fem::double0;
  double dpmd = fem::double0;
  double dpmm = fem::double0;
  float cad = fem::float0;
  int iim = fem::int0;
  int id1 = fem::int0;
  int id2 = fem::int0;
  float chi = fem::float0;
  double dbex = fem::double0;
  double dbey = fem::double0;
  double dbez = fem::double0;
  //C
  //C...Purpose: to generate timelike parton showers from given partons.
  //C
  //C...Initialization of cutoff masses etc.
  if (mstj(41) <= 0 || (mstj(41) == 1 && qmax <= parj(82)) ||
      qmax <= fem::min(parj(82), parj(83)) || mstj(41) >= 3) {
    return;
  }
  pmth(1, 21) = ulmass(cmn, 21);
  pmth(2, 21) = fem::sqrt(fem::pow2(pmth(1, 21)) + 0.25f * fem::pow2(parj(82)));
  pmth(3, 21) = 2.f * pmth(2, 21);
  pmth(4, 21) = pmth(3, 21);
  pmth(5, 21) = pmth(3, 21);
  pmth(1, 22) = ulmass(cmn, 22);
  pmth(2, 22) = fem::sqrt(fem::pow2(pmth(1, 22)) + 0.25f * fem::pow2(parj(83)));
  pmth(3, 22) = 2.f * pmth(2, 22);
  pmth(4, 22) = pmth(3, 22);
  pmth(5, 22) = pmth(3, 22);
  pmqth1 = parj(82);
  if (mstj(41) == 2) {
    pmqth1 = fem::min(parj(82), parj(83));
  }
  pmqth2 = pmth(2, 21);
  if (mstj(41) == 2) {
    pmqth2 = fem::min(pmth(2, 21), pmth(2, 22));
  }
  FEM_DO_SAFE(identifier_if, 1, 8) {
    pmth(1, identifier_if) = ulmass(cmn, identifier_if);
    pmth(2, identifier_if) = fem::sqrt(fem::pow2(pmth(1,
      identifier_if)) + 0.25f * fem::pow2(pmqth1));
    pmth(3, identifier_if) = pmth(2, identifier_if) + pmqth2;
    pmth(4, identifier_if) = fem::sqrt(fem::pow2(pmth(1,
      identifier_if)) + 0.25f * fem::pow2(parj(82))) + pmth(2, 21);
    pmth(5, identifier_if) = fem::sqrt(fem::pow2(pmth(1,
      identifier_if)) + 0.25f * fem::pow2(parj(83))) + pmth(2, 22);
  }
  pt2min = fem::pow2(fem::max(0.5f * parj(82), 1.1f * parj(81)));
  alams = fem::pow2(parj(81));
  alfm = fem::log(pt2min / alams);
  //C
  //C...Store positions of shower initiating partons.
  m3jc = 0;
  if (ip1 > 0 && ip1 <= fem::min(n, mstu(4) - mstu(32)) && ip2 == 0) {
    npa = 1;
    ipa(1) = ip1;
  }
  else if (fem::min(ip1, ip2) > 0 && fem::max(ip1, ip2) <= fem::min(n,
    mstu(4) - mstu(32))) {
    npa = 2;
    ipa(1) = ip1;
    ipa(2) = ip2;
  }
  else if (ip1 > 0 && ip1 <= fem::min(n, mstu(4) - mstu(32)) &&
    ip2 < 0 && ip2 >=  - 3) {
    npa = fem::iabs(ip2);
    FEM_DO_SAFE(i, 1, npa) {
      ipa(i) = ip1 + i - 1;
    }
  }
  else {
    luerrm(cmn, 12, "(LUSHOW:) failed to reconstruct showering system");
    if (mstu(21) >= 1) {
      return;
    }
  }
  //C
  //C...Check on phase space available for emission.
  irej = 0;
  FEM_DO_SAFE(j, 1, 5) {
    ps(j) = 0.f;
  }
  pm = 0.f;
  FEM_DO_SAFE(i, 1, npa) {
    kfla(i) = fem::iabs(k(ipa(i), 2));
    pma(i) = p(ipa(i), 5);
    if (kfla(i) != 0 && (kfla(i) <= 8 || kfla(i) == 21)) {
      pma(i) = pmth(3, kfla(i));
    }
    pm += pma(i);
    if (kfla(i) == 0 || (kfla(i) > 8 && kfla(i) != 21) || pma(i) > qmax) {
      irej++;
    }
    FEM_DO_SAFE(j, 1, 4) {
      ps(j) += p(ipa(i), j);
    }
  }
  if (irej == npa) {
    return;
  }
  ps(5) = fem::sqrt(fem::max(0.f, fem::pow2(ps(4)) - fem::pow2(ps(
    1)) - fem::pow2(ps(2)) - fem::pow2(ps(3))));
  if (npa == 1) {
    ps(5) = ps(4);
  }
  if (ps(5) <= pm + pmqth1) {
    return;
  }
  if (npa == 2 && mstj(47) >= 1) {
    if (kfla(1) >= 1 && kfla(1) <= 8 && kfla(2) >= 1 && kfla(2) <= 8) {
      m3jc = 1;
    }
    if (mstj(47) >= 2) {
      m3jc = 1;
    }
  }
  //C
  //C...Define imagined single initiator of shower for parton system.
  ns = n;
  if (n > mstu(4) - mstu(32) - 5) {
    luerrm(cmn, 11, "(LUSHOW:) no more memory left in LUJETS");
    if (mstu(21) >= 1) {
      return;
    }
  }
  if (npa >= 2) {
    k(n + 1, 1) = 11;
    k(n + 1, 2) = 21;
    k(n + 1, 3) = 0;
    k(n + 1, 4) = 0;
    k(n + 1, 5) = 0;
    p(n + 1, 1) = 0.f;
    p(n + 1, 2) = 0.f;
    p(n + 1, 3) = 0.f;
    p(n + 1, 4) = ps(5);
    p(n + 1, 5) = ps(5);
    v(n + 1, 5) = fem::pow2(ps(5));
    n++;
  }
  //C
  //C...Loop over partons that may branch.
  nep = npa;
  im = ns;
  if (npa == 1) {
    im = ns - 1;
  }
  statement_140:
  im++;
  if (n > ns) {
    if (im > n) {
      goto statement_380;
    }
    kflm = fem::iabs(k(im, 2));
    if (kflm == 0 || (kflm > 8 && kflm != 21)) {
      goto statement_140;
    }
    if (p(im, 5) < pmth(2, kflm)) {
      goto statement_140;
    }
    igm = k(im, 3);
  }
  else {
    igm = -1;
  }
  if (n + nep > mstu(4) - mstu(32) - 5) {
    luerrm(cmn, 11, "(LUSHOW:) no more memory left in LUJETS");
    if (mstu(21) >= 1) {
      return;
    }
  }
  //C
  //C...Position of aunt (sister to branching parton).
  //C...Origin and flavour of daughters.
  iau = 0;
  if (igm > 0) {
    if (k(im - 1, 3) == igm) {
      iau = im - 1;
    }
    if (n >= im + 1 && k(im + 1, 3) == igm) {
      iau = im + 1;
    }
  }
  if (igm >= 0) {
    k(im, 4) = n + 1;
    FEM_DO_SAFE(i, 1, nep) {
      k(n + i, 3) = im;
    }
  }
  else {
    k(n + 1, 3) = ipa(1);
  }
  if (igm <= 0) {
    FEM_DO_SAFE(i, 1, nep) {
      k(n + i, 2) = k(ipa(i), 2);
    }
  }
  else if (kflm != 21) {
    k(n + 1, 2) = k(im, 2);
    k(n + 2, 2) = k(im, 5);
  }
  else if (k(im, 5) == 21) {
    k(n + 1, 2) = 21;
    k(n + 2, 2) = 21;
  }
  else {
    k(n + 1, 2) = k(im, 5);
    k(n + 2, 2) = -k(im, 5);
  }
  //C
  //C...Reset flags on daughers and tries made.
  FEM_DO_SAFE(ip, 1, nep) {
    k(n + ip, 1) = 3;
    k(n + ip, 4) = 0;
    k(n + ip, 5) = 0;
    kfld(ip) = fem::iabs(k(n + ip, 2));
    itry(ip) = 0;
    isl(ip) = 0;
    isi(ip) = 0;
    if (kfld(ip) > 0 && (kfld(ip) <= 8 || kfld(ip) == 21)) {
      isi(ip) = 1;
    }
  }
  islm = 0;
  //C
  //C...Maximum virtuality of daughters.
  if (igm <= 0) {
    FEM_DO_SAFE(i, 1, npa) {
      if (npa >= 3) {
        p(n + i, 4) = (ps(4) * p(ipa(i), 4) - ps(1) * p(ipa(i), 1) -
          ps(2) * p(ipa(i), 2) - ps(3) * p(ipa(i), 3)) / ps(5);
      }
      p(n + i, 5) = fem::min(qmax, ps(5));
      if (npa >= 3) {
        p(n + i, 5) = fem::min(p(n + i, 5), p(n + i, 4));
      }
      if (isi(i) == 0) {
        p(n + i, 5) = p(ipa(i), 5);
      }
    }
  }
  else {
    if (mstj(43) <= 2) {
      pem = v(im, 2);
    }
    if (mstj(43) >= 3) {
      pem = p(im, 4);
    }
    p(n + 1, 5) = fem::min(p(im, 5), v(im, 1) * pem);
    p(n + 2, 5) = fem::min(p(im, 5), (1.f - v(im, 1)) * pem);
    if (k(n + 2, 2) == 22) {
      p(n + 2, 5) = pmth(1, 22);
    }
  }
  FEM_DO_SAFE(i, 1, nep) {
    pmsd(i) = p(n + i, 5);
    if (isi(i) == 1) {
      if (p(n + i, 5) <= pmth(3, kfld(i))) {
        p(n + i, 5) = pmth(1, kfld(i));
      }
    }
    v(n + i, 5) = fem::pow2(p(n + i, 5));
  }
  //C
  //C...Choose one of the daughters for evolution.
  statement_200:
  inum = 0;
  if (nep == 1) {
    inum = 1;
  }
  FEM_DO_SAFE(i, 1, nep) {
    if (inum == 0 && isl(i) == 1) {
      inum = i;
    }
  }
  FEM_DO_SAFE(i, 1, nep) {
    if (inum == 0 && itry(i) == 0 && isi(i) == 1) {
      if (p(n + i, 5) >= pmth(2, kfld(i))) {
        inum = i;
      }
    }
  }
  if (inum == 0) {
    rmax = 0.f;
    FEM_DO_SAFE(i, 1, nep) {
      if (isi(i) == 1 && pmsd(i) >= pmqth2) {
        rpm = p(n + i, 5) / pmsd(i);
        if (rpm > rmax && p(n + i, 5) >= pmth(2, kfld(i))) {
          rmax = rpm;
          inum = i;
        }
      }
    }
  }
  //C
  //C...Store information on choice of evolving daughter.
  inum = fem::max(1, inum);
  iep(1) = n + inum;
  FEM_DO_SAFE(i, 2, nep) {
    iep(i) = iep(i - 1) + 1;
    if (iep(i) > n + nep) {
      iep(i) = n + 1;
    }
  }
  FEM_DO_SAFE(i, 1, nep) {
    kfl(i) = fem::iabs(k(iep(i), 2));
  }
  itry(inum)++;
  if (itry(inum) > 200) {
    luerrm(cmn, 14, "(LUSHOW:) caught in infinite loop");
    if (mstu(21) >= 1) {
      return;
    }
  }
  z = 0.5f;
  if (kfl(1) == 0 || (kfl(1) > 8 && kfl(1) != 21)) {
    goto statement_300;
  }
  if (p(iep(1), 5) < pmth(2, kfl(1))) {
    goto statement_300;
  }
  //C
  //C...Calculate allowed z range.
  if (nep == 1) {
    pmed = ps(4);
  }
  else if (igm == 0 || mstj(43) <= 2) {
    pmed = p(im, 5);
  }
  else {
    if (inum == 1) {
      pmed = v(im, 1) * pem;
    }
    if (inum == 2) {
      pmed = (1.f - v(im, 1)) * pem;
    }
  }
  if (fem::mod(mstj(43), 2) == 1) {
    zc = pmth(2, 21) / pmed;
    zce = pmth(2, 22) / pmed;
  }
  else {
    zc = 0.5f * (1.f - fem::sqrt(fem::max(0.f, 1.f - fem::pow2((2.f * pmth(2,
      21) / pmed)))));
    if (zc < 1e-4f) {
      zc = fem::pow2((pmth(2, 21) / pmed));
    }
    zce = 0.5f * (1.f - fem::sqrt(fem::max(0.f, 1.f - fem::pow2((2.f * pmth(2,
      22) / pmed)))));
    if (zce < 1e-4f) {
      zce = fem::pow2((pmth(2, 22) / pmed));
    }
  }
  zc = fem::min(zc, 0.491f);
  zce = fem::min(zce, 0.491f);
  if ((mstj(41) == 1 && zc > 0.49f) || (mstj(41) == 2 && fem::min(zc,
      zce) > 0.49f)) {
    p(iep(1), 5) = pmth(1, kfl(1));
    v(iep(1), 5) = fem::pow2(p(iep(1), 5));
    goto statement_300;
  }
  //C
  //C...Integral of Altarelli-Parisi z kernel for QCD.
  if (mstj(49) == 0 && kfl(1) == 21) {
    fbr = 6.f * fem::log((1.f - zc) / zc) + mstj(45) * (0.5f - zc);
  }
  else if (mstj(49) == 0) {
    fbr = (8.f / 3.f) * fem::log((1.f - zc) / zc);
    //C
    //C...Integral of Altarelli-Parisi z kernel for scalar gluon.
  }
  else if (mstj(49) == 1 && kfl(1) == 21) {
    fbr = (parj(87) + mstj(45) * parj(88)) * (1.f - 2.f * zc);
  }
  else if (mstj(49) == 1) {
    fbr = (1.f - 2.f * zc) / 3.f;
    if (igm == 0 && m3jc == 1) {
      fbr = 4.f * fbr;
    }
    //C
    //C...Integral of Altarelli-Parisi z kernel for Abelian vector gluon.
  }
  else if (kfl(1) == 21) {
    fbr = 6.f * mstj(45) * (0.5f - zc);
  }
  else {
    fbr = 2.f * fem::log((1.f - zc) / zc);
  }
  //C
  //C...Integral of Altarelli-Parisi kernel for photon emission.
  if (mstj(41) == 2 && kfl(1) >= 1 && kfl(1) <= 8) {
    fbre = fem::pow2((kchg(kfl(1), 1) / 3.f)) * 2.f * fem::log((1.f -
      zce) / zce);
  }
  //C
  //C...Inner veto algorithm starts. Find maximum mass for evolution.
  statement_260:
  pms = v(iep(1), 5);
  if (igm >= 0) {
    pm2 = 0.f;
    FEM_DO_SAFE(i, 2, nep) {
      pm = p(iep(i), 5);
      if (kfl(i) > 0 && (kfl(i) <= 8 || kfl(i) == 21)) {
        pm = pmth(2, kfl(i));
      }
      pm2 += pm;
    }
    pms = fem::min(pms, fem::pow2((p(im, 5) - pm2)));
  }
  //C
  //C...Select mass for daughter in QCD evolution.
  b0 = 27.f / 6.f;
  FEM_DO_SAFE(identifier_if, 4, mstj(45)) {
    if (pms > 4.f * fem::pow2(pmth(2, identifier_if))) {
      b0 = (33.f - 2.f * identifier_if) / 6.f;
    }
  }
  if (mstj(44) <= 0) {
    pmsqcd = pms * fem::exp(fem::max(-100.f, fem::log(rlu(cmn, 0)) *
      paru(2) / (paru(111) * fbr)));
  }
  else if (mstj(44) == 1) {
    pmsqcd = 4.f * alams * fem::pow((0.25f * pms / alams), (fem::pow(rlu(cmn,
      0), (b0 / fbr))));
  }
  else {
    pmsqcd = pms * fem::pow(rlu(cmn, 0), (alfm * b0 / fbr));
  }
  if (zc > 0.49f || pmsqcd <= fem::pow2(pmth(4, kfl(1)))) {
    pmsqcd = fem::pow2(pmth(2, kfl(1)));
  }
  v(iep(1), 5) = pmsqcd;
  mce = 1;
  //C
  //C...Select mass for daughter in QED evolution.
  if (mstj(41) == 2 && kfl(1) >= 1 && kfl(1) <= 8) {
    pmsqed = pms * fem::exp(fem::max(-100.f, fem::log(rlu(cmn, 0)) *
      paru(2) / (paru(101) * fbre)));
    if (zce > 0.49f || pmsqed <= fem::pow2(pmth(5, kfl(1)))) {
      pmsqed = fem::pow2(pmth(2, kfl(1)));
    }
    if (pmsqed > pmsqcd) {
      v(iep(1), 5) = pmsqed;
      mce = 2;
    }
  }
  //C
  //C...Check whether daughter mass below cutoff.
  p(iep(1), 5) = fem::sqrt(v(iep(1), 5));
  if (p(iep(1), 5) <= pmth(3, kfl(1))) {
    p(iep(1), 5) = pmth(1, kfl(1));
    v(iep(1), 5) = fem::pow2(p(iep(1), 5));
    goto statement_300;
  }
  //C
  //C...Select z value of branching: q -> qgamma.
  if (mce == 2) {
    z = 1.f - (1.f - zce) * fem::pow((zce / (1.f - zce)), rlu(cmn, 0));
    if (1.f + fem::pow2(z) < 2.f * rlu(cmn, 0)) {
      goto statement_260;
    }
    k(iep(1), 5) = 22;
    //C
    //C...Select z value of branching: q -> qg, g -> gg, g -> qqbar.
  }
  else if (mstj(49) != 1 && kfl(1) != 21) {
    z = 1.f - (1.f - zc) * fem::pow((zc / (1.f - zc)), rlu(cmn, 0));
    if (1.f + fem::pow2(z) < 2.f * rlu(cmn, 0)) {
      goto statement_260;
    }
    k(iep(1), 5) = 21;
  }
  else if (mstj(49) == 0 && mstj(45) * (0.5f - zc) < rlu(cmn, 0) * fbr) {
    z = (1.f - zc) * fem::pow((zc / (1.f - zc)), rlu(cmn, 0));
    if (rlu(cmn, 0) > 0.5f) {
      z = 1.f - z;
    }
    if (fem::pow2((1.f - z * (1.f - z))) < rlu(cmn, 0)) {
      goto statement_260;
    }
    k(iep(1), 5) = 21;
  }
  else if (mstj(49) != 1) {
    z = zc + (1.f - 2.f * zc) * rlu(cmn, 0);
    if (fem::pow2(z) + fem::pow2((1.f - z)) < rlu(cmn, 0)) {
      goto statement_260;
    }
    kflb = 1 + fem::fint(mstj(45) * rlu(cmn, 0));
    pmq = 4.f * fem::pow2(pmth(2, kflb)) / v(iep(1), 5);
    if (pmq >= 1.f) {
      goto statement_260;
    }
    pmq0 = 4.f * fem::pow2(pmth(2, 21)) / v(iep(1), 5);
    if (fem::mod(mstj(43), 2) == 0 && (1.f + 0.5f * pmq) * fem::sqrt(
        1.f - pmq) < rlu(cmn, 0) * (1.f + 0.5f * pmq0) * fem::sqrt(
        1.f - pmq0)) {
      goto statement_260;
    }
    k(iep(1), 5) = kflb;
    //C
    //C...Ditto for scalar gluon model.
  }
  else if (kfl(1) != 21) {
    z = 1.f - fem::sqrt(fem::pow2(zc) + rlu(cmn, 0) * (1.f - 2.f * zc));
    k(iep(1), 5) = 21;
  }
  else if (rlu(cmn, 0) * (parj(87) + mstj(45) * parj(88)) <= parj(87)) {
    z = zc + (1.f - 2.f * zc) * rlu(cmn, 0);
    k(iep(1), 5) = 21;
  }
  else {
    z = zc + (1.f - 2.f * zc) * rlu(cmn, 0);
    kflb = 1 + fem::fint(mstj(45) * rlu(cmn, 0));
    pmq = 4.f * fem::pow2(pmth(2, kflb)) / v(iep(1), 5);
    if (pmq >= 1.f) {
      goto statement_260;
    }
    k(iep(1), 5) = kflb;
  }
  if (mce == 1 && mstj(44) >= 2) {
    if (z * (1.f - z) * v(iep(1), 5) < pt2min) {
      goto statement_260;
    }
    if (alfm / fem::log(v(iep(1), 5) * z * (1.f - z) / alams) < rlu(cmn, 0)) {
      goto statement_260;
    }
  }
  //C
  //C...Check if z consistent with chosen m.
  if (kfl(1) == 21) {
    kflgd1 = fem::iabs(k(iep(1), 5));
    kflgd2 = kflgd1;
  }
  else {
    kflgd1 = kfl(1);
    kflgd2 = fem::iabs(k(iep(1), 5));
  }
  if (nep == 1) {
    ped = ps(4);
  }
  else if (nep >= 3) {
    ped = p(iep(1), 4);
  }
  else if (igm == 0 || mstj(43) <= 2) {
    ped = 0.5f * (v(im, 5) + v(iep(1), 5) - fem::pow2(pm2)) / p(im, 5);
  }
  else {
    if (iep(1) == n + 1) {
      ped = v(im, 1) * pem;
    }
    if (iep(1) == n + 2) {
      ped = (1.f - v(im, 1)) * pem;
    }
  }
  if (fem::mod(mstj(43), 2) == 1) {
    pmqth3 = 0.5f * parj(82);
    if (kflgd2 == 22) {
      pmqth3 = 0.5f * parj(83);
    }
    pmq1 = (fem::pow2(pmth(1, kflgd1)) + fem::pow2(pmqth3)) / v(iep(1), 5);
    pmq2 = (fem::pow2(pmth(1, kflgd2)) + fem::pow2(pmqth3)) / v(iep(1), 5);
    zd = fem::sqrt(fem::max(0.f, (1.f - v(iep(1), 5) / fem::pow2(
      ped)) * (fem::pow2((1.f - pmq1 - pmq2)) - 4.f * pmq1 * pmq2)));
    zh = 1.f + pmq1 - pmq2;
  }
  else {
    zd = fem::sqrt(fem::max(0.f, 1.f - v(iep(1), 5) / fem::pow2(ped)));
    zh = 1.f;
  }
  zl = 0.5f * (zh - zd);
  zu = 0.5f * (zh + zd);
  if (z < zl || z > zu) {
    goto statement_260;
  }
  if (kfl(1) == 21) {
    v(iep(1), 3) = fem::log(zu * (1.f - zl) / fem::max(1e-20f, zl * (
      1.f - zu)));
  }
  if (kfl(1) != 21) {
    v(iep(1), 3) = fem::log((1.f - zl) / fem::max(1e-10f, 1.f - zu));
  }
  //C
  //C...Three-jet matrix element correction.
  if (igm == 0 && m3jc == 1) {
    x1 = z * (1.f + v(iep(1), 5) / v(ns + 1, 5));
    x2 = 1.f - v(iep(1), 5) / v(ns + 1, 5);
    x3 = (1.f - x1) + (1.f - x2);
    if (mce == 2) {
      ki1 = k(ipa(inum), 2);
      ki2 = k(ipa(3 - inum), 2);
      qf1 = kchg(fem::iabs(ki1), 1) * fem::isign(1, ki1) / 3.f;
      qf2 = kchg(fem::iabs(ki2), 1) * fem::isign(1, ki2) / 3.f;
      wshow = fem::pow2(qf1) * (1.f - x1) / x3 * (1.f + fem::pow2((
        x1 / (2.f - x2)))) + fem::pow2(qf2) * (1.f - x2) / x3 * (
        1.f + fem::pow2((x2 / (2.f - x1))));
      wme = fem::pow2((qf1 * (1.f - x1) / x3 - qf2 * (1.f - x2) /
        x3)) * (fem::pow2(x1) + fem::pow2(x2));
    }
    else if (mstj(49) != 1) {
      wshow = 1.f + (1.f - x1) / x3 * fem::pow2((x1 / (2.f - x2))) + (
        1.f - x2) / x3 * fem::pow2((x2 / (2.f - x1)));
      wme = fem::pow2(x1) + fem::pow2(x2);
    }
    else {
      wshow = 4.f * x3 * ((1.f - x1) / fem::pow2((2.f - x2)) + (1.f -
        x2) / fem::pow2((2.f - x1)));
      wme = fem::pow2(x3);
    }
    if (wme < rlu(cmn, 0) * wshow) {
      goto statement_260;
    }
    //C
    //C...Impose angular ordering by rejection of nonordered emission.
  }
  else if (mce == 1 && igm > 0 && mstj(42) >= 2) {
    maom = 1;
    zm = v(im, 1);
    if (iep(1) == n + 2) {
      zm = 1.f - v(im, 1);
    }
    the2id = z * (1.f - z) * fem::pow2((zm * p(im, 4))) / v(iep(1), 5);
    iaom = im;
    statement_290:
    if (k(iaom, 5) == 22) {
      iaom = k(iaom, 3);
      if (k(iaom, 3) <= ns) {
        maom = 0;
      }
      if (maom == 1) {
        goto statement_290;
      }
    }
    if (maom == 1) {
      the2im = v(iaom, 1) * (1.f - v(iaom, 1)) * fem::pow2(p(iaom,
        4)) / v(iaom, 5);
      if (the2id < the2im) {
        goto statement_260;
      }
    }
  }
  //C
  //C...Impose user-defined maximum angle at first branching.
  if (mstj(48) == 1) {
    if (nep == 1 && im == ns) {
      the2id = z * (1.f - z) * fem::pow2(ps(4)) / v(iep(1), 5);
      if (the2id < 1.f / fem::pow2(parj(85))) {
        goto statement_260;
      }
    }
    else if (nep == 2 && iep(1) == ns + 2) {
      the2id = z * (1.f - z) * fem::pow2((0.5f * p(im, 4))) / v(iep(1), 5);
      if (the2id < 1.f / fem::pow2(parj(85))) {
        goto statement_260;
      }
    }
    else if (nep == 2 && iep(1) == ns + 3) {
      the2id = z * (1.f - z) * fem::pow2((0.5f * p(im, 4))) / v(iep(1), 5);
      if (the2id < 1.f / fem::pow2(parj(86))) {
        goto statement_260;
      }
    }
  }
  //C
  //C...End of inner veto algorithm. Check if only one leg evolved so far.
  statement_300:
  v(iep(1), 1) = z;
  isl(1) = 0;
  isl(2) = 0;
  if (nep == 1) {
    goto statement_330;
  }
  if (nep == 2 && p(iep(1), 5) + p(iep(2), 5) >= p(im, 5)) {
    goto statement_200;
  }
  FEM_DO_SAFE(i, 1, nep) {
    if (itry(i) == 0 && kfld(i) > 0 && (kfld(i) <= 8 || kfld(i) == 21)) {
      if (p(n + i, 5) >= pmth(2, kfld(i))) {
        goto statement_200;
      }
    }
  }
  //C
  //C...Check if chosen multiplet m1,m2,z1,z2 is physical.
  if (nep == 3) {
    pa1s = (p(n + 1, 4) + p(n + 1, 5)) * (p(n + 1, 4) - p(n + 1, 5));
    pa2s = (p(n + 2, 4) + p(n + 2, 5)) * (p(n + 2, 4) - p(n + 2, 5));
    pa3s = (p(n + 3, 4) + p(n + 3, 5)) * (p(n + 3, 4) - p(n + 3, 5));
    pts = 0.25f * (2.f * pa1s * pa2s + 2.f * pa1s * pa3s + 2.f *
      pa2s * pa3s - fem::pow2(pa1s) - fem::pow2(pa2s) - fem::pow2(
      pa3s)) / pa1s;
    if (pts <= 0.f) {
      goto statement_200;
    }
  }
  else if (igm == 0 || mstj(43) <= 2 || fem::mod(mstj(43), 2) == 0) {
    FEM_DO_SAFE(i1, n + 1, n + 2) {
      kflda = fem::iabs(k(i1, 2));
      if (kflda == 0 || (kflda > 8 && kflda != 21)) {
        goto statement_320;
      }
      if (p(i1, 5) < pmth(2, kflda)) {
        goto statement_320;
      }
      if (kflda == 21) {
        kflgd1 = fem::iabs(k(i1, 5));
        kflgd2 = kflgd1;
      }
      else {
        kflgd1 = kflda;
        kflgd2 = fem::iabs(k(i1, 5));
      }
      i2 = 2 * n + 3 - i1;
      if (igm == 0 || mstj(43) <= 2) {
        ped = 0.5f * (v(im, 5) + v(i1, 5) - v(i2, 5)) / p(im, 5);
      }
      else {
        if (i1 == n + 1) {
          zm = v(im, 1);
        }
        if (i1 == n + 2) {
          zm = 1.f - v(im, 1);
        }
        pml = fem::sqrt(fem::pow2((v(im, 5) - v(n + 1, 5) - v(n + 2,
          5))) - 4.f * v(n + 1, 5) * v(n + 2, 5));
        ped = pem * (0.5f * (v(im, 5) - pml + v(i1, 5) - v(i2, 5)) +
          pml * zm) / v(im, 5);
      }
      if (fem::mod(mstj(43), 2) == 1) {
        pmqth3 = 0.5f * parj(82);
        if (kflgd2 == 22) {
          pmqth3 = 0.5f * parj(83);
        }
        pmq1 = (fem::pow2(pmth(1, kflgd1)) + fem::pow2(pmqth3)) / v(i1, 5);
        pmq2 = (fem::pow2(pmth(1, kflgd2)) + fem::pow2(pmqth3)) / v(i1, 5);
        zd = fem::sqrt(fem::max(0.f, (1.f - v(i1, 5) / fem::pow2(
          ped)) * (fem::pow2((1.f - pmq1 - pmq2)) - 4.f * pmq1 *
          pmq2)));
        zh = 1.f + pmq1 - pmq2;
      }
      else {
        zd = fem::sqrt(fem::max(0.f, 1.f - v(i1, 5) / fem::pow2(ped)));
        zh = 1.f;
      }
      zl = 0.5f * (zh - zd);
      zu = 0.5f * (zh + zd);
      if (i1 == n + 1 && (v(i1, 1) < zl || v(i1, 1) > zu)) {
        isl(1) = 1;
      }
      if (i1 == n + 2 && (v(i1, 1) < zl || v(i1, 1) > zu)) {
        isl(2) = 1;
      }
      if (kflda == 21) {
        v(i1, 4) = fem::log(zu * (1.f - zl) / fem::max(1e-20f, zl * (
          1.f - zu)));
      }
      if (kflda != 21) {
        v(i1, 4) = fem::log((1.f - zl) / fem::max(1e-10f, 1.f - zu));
      }
      statement_320:;
    }
    if (isl(1) == 1 && isl(2) == 1 && islm != 0) {
      isl(3 - islm) = 0;
      islm = 3 - islm;
    }
    else if (isl(1) == 1 && isl(2) == 1) {
      zdr1 = fem::max(0.f, v(n + 1, 3) / v(n + 1, 4) - 1.f);
      zdr2 = fem::max(0.f, v(n + 2, 3) / v(n + 2, 4) - 1.f);
      if (zdr2 > rlu(cmn, 0) * (zdr1 + zdr2)) {
        isl(1) = 0;
      }
      if (isl(1) == 1) {
        isl(2) = 0;
      }
      if (isl(1) == 0) {
        islm = 1;
      }
      if (isl(2) == 0) {
        islm = 2;
      }
    }
    if (isl(1) == 1 || isl(2) == 1) {
      goto statement_200;
    }
  }
  if (igm > 0 && fem::mod(mstj(43), 2) == 1 && (p(n + 1, 5) >= pmth(2,
      kfld(1)) || p(n + 2, 5) >= pmth(2, kfld(2)))) {
    pmq1 = v(n + 1, 5) / v(im, 5);
    pmq2 = v(n + 2, 5) / v(im, 5);
    zd = fem::sqrt(fem::max(0.f, (1.f - v(im, 5) / fem::pow2(pem)) * (
      fem::pow2((1.f - pmq1 - pmq2)) - 4.f * pmq1 * pmq2)));
    zh = 1.f + pmq1 - pmq2;
    zl = 0.5f * (zh - zd);
    zu = 0.5f * (zh + zd);
    if (v(im, 1) < zl || v(im, 1) > zu) {
      goto statement_200;
    }
  }
  //C
  //C...Accepted branch. Construct four-momentum for initial partons.
  statement_330:
  mazip = 0;
  mazic = 0;
  if (nep == 1) {
    p(n + 1, 1) = 0.f;
    p(n + 1, 2) = 0.f;
    p(n + 1, 3) = fem::sqrt(fem::max(0.f, (p(ipa(1), 4) + p(n + 1,
      5)) * (p(ipa(1), 4) - p(n + 1, 5))));
    p(n + 1, 4) = p(ipa(1), 4);
    v(n + 1, 2) = p(n + 1, 4);
  }
  else if (igm == 0 && nep == 2) {
    ped1 = 0.5f * (v(im, 5) + v(n + 1, 5) - v(n + 2, 5)) / p(im, 5);
    p(n + 1, 1) = 0.f;
    p(n + 1, 2) = 0.f;
    p(n + 1, 3) = fem::sqrt(fem::max(0.f, (ped1 + p(n + 1, 5)) * (
      ped1 - p(n + 1, 5))));
    p(n + 1, 4) = ped1;
    p(n + 2, 1) = 0.f;
    p(n + 2, 2) = 0.f;
    p(n + 2, 3) = -p(n + 1, 3);
    p(n + 2, 4) = p(im, 5) - ped1;
    v(n + 1, 2) = p(n + 1, 4);
    v(n + 2, 2) = p(n + 2, 4);
  }
  else if (nep == 3) {
    p(n + 1, 1) = 0.f;
    p(n + 1, 2) = 0.f;
    p(n + 1, 3) = fem::sqrt(fem::max(0.f, pa1s));
    p(n + 2, 1) = fem::sqrt(pts);
    p(n + 2, 2) = 0.f;
    p(n + 2, 3) = 0.5f * (pa3s - pa2s - pa1s) / p(n + 1, 3);
    p(n + 3, 1) = -p(n + 2, 1);
    p(n + 3, 2) = 0.f;
    p(n + 3, 3) = -(p(n + 1, 3) + p(n + 2, 3));
    v(n + 1, 2) = p(n + 1, 4);
    v(n + 2, 2) = p(n + 2, 4);
    v(n + 3, 2) = p(n + 3, 4);
    //C
    //C...Construct transverse momentum for ordinary branching in shower.
  }
  else {
    zm = v(im, 1);
    pzm = fem::sqrt(fem::max(0.f, (pem + p(im, 5)) * (pem - p(im, 5))));
    pmls = fem::pow2((v(im, 5) - v(n + 1, 5) - v(n + 2, 5))) - 4.f * v(n + 1,
      5) * v(n + 2, 5);
    if (pzm <= 0.f) {
      pts = 0.f;
    }
    else if (fem::mod(mstj(43), 2) == 1) {
      pts = (fem::pow2(pem) * (zm * (1.f - zm) * v(im, 5) - (1.f -
        zm) * v(n + 1, 5) - zm * v(n + 2, 5)) - 0.25f * pmls) /
        fem::pow2(pzm);
    }
    else {
      pts = pmls * (zm * (1.f - zm) * fem::pow2(pem) / v(im, 5) -
        0.25f) / fem::pow2(pzm);
    }
    pt = fem::sqrt(fem::max(0.f, pts));
    //C
    //C...Find coefficient of azimuthal asymmetry due to gluon polarization.
    hazip = 0.f;
    if (mstj(49) != 1 && fem::mod(mstj(46), 2) == 1 && k(im,
        2) == 21 && iau != 0) {
      if (k(igm, 3) != 0) {
        mazip = 1;
      }
      zau = v(igm, 1);
      if (iau == im + 1) {
        zau = 1.f - v(igm, 1);
      }
      if (mazip == 0) {
        zau = 0.f;
      }
      if (k(igm, 2) != 21) {
        hazip = 2.f * zau / (1.f + fem::pow2(zau));
      }
      else {
        hazip = fem::pow2((zau / (1.f - zau * (1.f - zau))));
      }
      if (k(n + 1, 2) != 21) {
        hazip = hazip * (-2.f * zm * (1.f - zm)) / (1.f - 2.f * zm * (
          1.f - zm));
      }
      else {
        hazip = hazip * fem::pow2((zm * (1.f - zm) / (1.f - zm * (1.f - zm))));
      }
    }
    //C
    //C...Find coefficient of azimuthal asymmetry due to soft gluon
    //C...interference.
    hazic = 0.f;
    if (mstj(46) >= 2 && (k(n + 1, 2) == 21 || k(n + 2, 2) == 21) && iau != 0) {
      if (k(igm, 3) != 0) {
        mazic = n + 1;
      }
      if (k(igm, 3) != 0 && k(n + 1, 2) != 21) {
        mazic = n + 2;
      }
      if (k(igm, 3) != 0 && k(n + 1, 2) == 21 && k(n + 2, 2) == 21 &&
          zm > 0.5f) {
        mazic = n + 2;
      }
      if (k(iau, 2) == 22) {
        mazic = 0;
      }
      zs = zm;
      if (mazic == n + 2) {
        zs = 1.f - zm;
      }
      zgm = v(igm, 1);
      if (iau == im - 1) {
        zgm = 1.f - v(igm, 1);
      }
      if (mazic == 0) {
        zgm = 1.f;
      }
      hazic = (p(im, 5) / p(igm, 5)) * fem::sqrt((1.f - zs) * (1.f -
        zgm) / (zs * zgm));
      hazic = fem::min(0.95f, hazic);
    }
  }
  //C
  //C...Construct kinematics for ordinary branching in shower.
  statement_340:
  if (nep == 2 && igm > 0) {
    if (fem::mod(mstj(43), 2) == 1) {
      p(n + 1, 4) = pem * v(im, 1);
    }
    else {
      p(n + 1, 4) = pem * (0.5f * (v(im, 5) - fem::sqrt(pmls) + v(n + 1,
        5) - v(n + 2, 5)) + fem::sqrt(pmls) * zm) / v(im, 5);
    }
    phi = paru(2) * rlu(cmn, 0);
    p(n + 1, 1) = pt * fem::cos(phi);
    p(n + 1, 2) = pt * fem::sin(phi);
    if (pzm > 0.f) {
      p(n + 1, 3) = 0.5f * (v(n + 2, 5) - v(n + 1, 5) - v(im, 5) +
        2.f * pem * p(n + 1, 4)) / pzm;
    }
    else {
      p(n + 1, 3) = 0.f;
    }
    p(n + 2, 1) = -p(n + 1, 1);
    p(n + 2, 2) = -p(n + 1, 2);
    p(n + 2, 3) = pzm - p(n + 1, 3);
    p(n + 2, 4) = pem - p(n + 1, 4);
    if (mstj(43) <= 2) {
      v(n + 1, 2) = (pem * p(n + 1, 4) - pzm * p(n + 1, 3)) / p(im, 5);
      v(n + 2, 2) = (pem * p(n + 2, 4) - pzm * p(n + 2, 3)) / p(im, 5);
    }
  }
  //C
  //C...Rotate and boost daughters.
  if (igm > 0) {
    if (mstj(43) <= 2) {
      bex = p(igm, 1) / p(igm, 4);
      bey = p(igm, 2) / p(igm, 4);
      bez = p(igm, 3) / p(igm, 4);
      ga = p(igm, 4) / p(igm, 5);
      gabep = ga * (ga * (bex * p(im, 1) + bey * p(im, 2) + bez * p(im,
        3)) / (1.f + ga) - p(im, 4));
    }
    else {
      bex = 0.f;
      bey = 0.f;
      bez = 0.f;
      ga = 1.f;
      gabep = 0.f;
    }
    the = ulangl(cmn, p(im, 3) + gabep * bez, fem::sqrt(fem::pow2((p(im,
      1) + gabep * bex)) + fem::pow2((p(im, 2) + gabep * bey))));
    phi = ulangl(cmn, p(im, 1) + gabep * bex, p(im, 2) + gabep * bey);
    FEM_DO_SAFE(i, n + 1, n + 2) {
      dp(1) = fem::dble(fem::cos(the) * fem::cos(phi) * p(i, 1) -
        fem::sin(phi) * p(i, 2) + fem::sin(the) * fem::cos(phi) * p(i,
        3));
      dp(2) = fem::dble(fem::cos(the) * fem::sin(phi) * p(i, 1) +
        fem::cos(phi) * p(i, 2) + fem::sin(the) * fem::sin(phi) * p(i,
        3));
      dp(3) = fem::dble(-fem::sin(the) * p(i, 1) + fem::cos(the) * p(i, 3));
      dp(4) = fem::dble(p(i, 4));
      dbp = fem::dble(bex) * dp(1) + fem::dble(bey) * dp(2) +
        fem::dble(bez) * dp(3);
      dgabp = fem::dble(ga) * (fem::dble(ga) * dbp / (1e0 + fem::dble(
        ga)) + dp(4));
      p(i, 1) = fem::sngl(dp(1) + dgabp * fem::dble(bex));
      p(i, 2) = fem::sngl(dp(2) + dgabp * fem::dble(bey));
      p(i, 3) = fem::sngl(dp(3) + dgabp * fem::dble(bez));
      p(i, 4) = ga * fem::sngl(dp(4) + dbp);
    }
  }
  //C
  //C...Weight with azimuthal distribution, if required.
  if (mazip != 0 || mazic != 0) {
    FEM_DO_SAFE(j, 1, 3) {
      dpt(1, j) = fem::dble(p(im, j));
      dpt(2, j) = fem::dble(p(iau, j));
      dpt(3, j) = fem::dble(p(n + 1, j));
    }
    dpma = dpt(1, 1) * dpt(2, 1) + dpt(1, 2) * dpt(2, 2) + dpt(1,
      3) * dpt(2, 3);
    dpmd = dpt(1, 1) * dpt(3, 1) + dpt(1, 2) * dpt(3, 2) + dpt(1,
      3) * dpt(3, 3);
    dpmm = fem::pow2(dpt(1, 1)) + fem::pow2(dpt(1, 2)) + fem::pow2(dpt(1, 3));
    FEM_DO_SAFE(j, 1, 3) {
      dpt(4, j) = dpt(2, j) - dpma * dpt(1, j) / dpmm;
      dpt(5, j) = dpt(3, j) - dpmd * dpt(1, j) / dpmm;
    }
    dpt(4, 4) = fem::dsqrt(fem::pow2(dpt(4, 1)) + fem::pow2(dpt(4,
      2)) + fem::pow2(dpt(4, 3)));
    dpt(5, 4) = fem::dsqrt(fem::pow2(dpt(5, 1)) + fem::pow2(dpt(5,
      2)) + fem::pow2(dpt(5, 3)));
    //Clin-5/2012:
    //C        IF(MIN(DPT(4,4),DPT(5,4)).GT.0.1*PARJ(82)) THEN
    if (fem::sngl(fem::min(dpt(4, 4), dpt(5, 4))) > (0.1f * parj(82))) {
      cad = fem::sngl((dpt(4, 1) * dpt(5, 1) + dpt(4, 2) * dpt(5,
        2) + dpt(4, 3) * dpt(5, 3)) / (dpt(4, 4) * dpt(5, 4)));
      if (mazip != 0) {
        if (1.f + hazip * (2.f * fem::pow2(cad) - 1.f) < rlu(cmn,
            0) * (1.f + fem::abs(hazip))) {
          goto statement_340;
        }
      }
      if (mazic != 0) {
        if (mazic == n + 2) {
          cad = -cad;
        }
        if ((1.f - hazic) * (1.f - hazic * cad) / (1.f + fem::pow2(
            hazic) - 2.f * hazic * cad) < rlu(cmn, 0)) {
          goto statement_340;
        }
      }
    }
  }
  //C
  //C...Continue loop over partons that may branch, until none left.
  if (igm >= 0) {
    k(im, 1) = 14;
  }
  n += nep;
  nep = 2;
  if (n > mstu(4) - mstu(32) - 5) {
    luerrm(cmn, 11, "(LUSHOW:) no more memory left in LUJETS");
    if (mstu(21) >= 1) {
      n = ns;
    }
    if (mstu(21) >= 1) {
      return;
    }
  }
  goto statement_140;
  //C
  //C...Set information on imagined shower initiator.
  statement_380:
  if (npa >= 2) {
    k(ns + 1, 1) = 11;
    k(ns + 1, 2) = 94;
    k(ns + 1, 3) = ip1;
    if (ip2 > 0 && ip2 < ip1) {
      k(ns + 1, 3) = ip2;
    }
    k(ns + 1, 4) = ns + 2;
    k(ns + 1, 5) = ns + 1 + npa;
    iim = 1;
  }
  else {
    iim = 0;
  }
  //C
  //C...Reconstruct string drawing information.
  FEM_DO_SAFE(i, ns + 1 + iim, n) {
    if (k(i, 1) <= 10 && k(i, 2) == 22) {
      k(i, 1) = 1;
    }
    else if (k(i, 1) <= 10) {
      k(i, 4) = mstu(5) * (k(i, 4) / mstu(5));
      k(i, 5) = mstu(5) * (k(i, 5) / mstu(5));
    }
    else if (k(fem::mod(k(i, 4), mstu(5)) + 1, 2) != 22) {
      id1 = fem::mod(k(i, 4), mstu(5));
      if (k(i, 2) >= 1 && k(i, 2) <= 8) {
        id1 = fem::mod(k(i, 4), mstu(5)) + 1;
      }
      id2 = 2 * fem::mod(k(i, 4), mstu(5)) + 1 - id1;
      k(i, 4) = mstu(5) * (k(i, 4) / mstu(5)) + id1;
      k(i, 5) = mstu(5) * (k(i, 5) / mstu(5)) + id2;
      k(id1, 4) += mstu(5) * i;
      k(id1, 5) += mstu(5) * id2;
      k(id2, 4) += mstu(5) * id1;
      k(id2, 5) += mstu(5) * i;
    }
    else {
      id1 = fem::mod(k(i, 4), mstu(5));
      id2 = id1 + 1;
      k(i, 4) = mstu(5) * (k(i, 4) / mstu(5)) + id1;
      k(i, 5) = mstu(5) * (k(i, 5) / mstu(5)) + id1;
      k(id1, 4) += mstu(5) * i;
      k(id1, 5) += mstu(5) * i;
      k(id2, 4) = 0;
      k(id2, 5) = 0;
    }
  }
  //C
  //C...Transformation from CM frame.
  if (npa >= 2) {
    bex = ps(1) / ps(4);
    bey = ps(2) / ps(4);
    bez = ps(3) / ps(4);
    ga = ps(4) / ps(5);
    gabep = ga * (ga * (bex * p(ipa(1), 1) + bey * p(ipa(1), 2) +
      bez * p(ipa(1), 3)) / (1.f + ga) - p(ipa(1), 4));
  }
  else {
    bex = 0.f;
    bey = 0.f;
    bez = 0.f;
    gabep = 0.f;
  }
  the = ulangl(cmn, p(ipa(1), 3) + gabep * bez, fem::sqrt(fem::pow2((p(ipa(1),
    1) + gabep * bex)) + fem::pow2((p(ipa(1), 2) + gabep * bey))));
  phi = ulangl(cmn, p(ipa(1), 1) + gabep * bex, p(ipa(1), 2) + gabep * bey);
  if (npa == 3) {
    chi = ulangl(cmn, fem::cos(the) * fem::cos(phi) * (p(ipa(2), 1) +
      gabep * bex) + fem::cos(the) * fem::sin(phi) * (p(ipa(2), 2) +
      gabep * bey) - fem::sin(the) * (p(ipa(2), 3) + gabep * bez),
      -fem::sin(phi) * (p(ipa(2), 1) + gabep * bex) + fem::cos(phi) * (
      p(ipa(2), 2) + gabep * bey));
    ludbrb(ns + 1, n, 0.f, chi, 0e0, 0e0, 0e0);
  }
  dbex = fem::dble(bex);
  dbey = fem::dble(bey);
  dbez = fem::dble(bez);
  ludbrb(ns + 1, n, the, phi, dbex, dbey, dbez);
  //C
  //C...Decay vertex of shower.
  FEM_DO_SAFE(i, ns + 1, n) {
    FEM_DO_SAFE(j, 1, 5) {
      v(i, j) = v(ip1, j);
    }
  }
  //C
  //C...Delete trivial shower, else connect initiators.
  if (n == ns + npa + iim) {
    n = ns;
  }
  else {
    FEM_DO_SAFE(ip, 1, npa) {
      k(ipa(ip), 1) = 14;
      k(ipa(ip), 4) += ns + iim + ip;
      k(ipa(ip), 5) += ns + iim + ip;
      k(ns + iim + ip, 3) = ipa(ip);
      if (iim == 1 && mstu(16) != 2) {
        k(ns + iim + ip, 3) = ns + 1;
      }
      k(ns + iim + ip, 4) += mstu(5) * ipa(ip);
      k(ns + iim + ip, 5) += mstu(5) * ipa(ip);
    }
  }
  //C
}

struct luboei_save
{
  arr<int> kfbe;

  luboei_save() :
    kfbe(dimension(9), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
luboei(
  common& cmn,
  int const& nsav)
{
  FEM_CMN_SVE(luboei);
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  //
  arr_ref<int> kfbe(sve.kfbe, dimension(9));
  if (is_called_first_time) {
    static const int values[] = {
      211, -211, 111, 321, -321, 130, 310, 221, 331
    };
    fem::data_of_type<int>(FEM_VALUES_AND_SIZE),
      kfbe;
  }
  int j = fem::int0;
  arr_1d<4, double> dps(fem::fill0);
  int i = fem::int0;
  float pecm = fem::float0;
  arr_1d<10, int> nbe(dim1(0, 9), fem::fill0);
  int ibe = fem::int0;
  float pmhq = fem::float0;
  float qdel = fem::float0;
  int nbin = fem::int0;
  float beex = fem::float0;
  float bert = fem::float0;
  int ibin = fem::int0;
  float qbin = fem::float0;
  arr_1d<100, float> bei(fem::fill0);
  int i1m = fem::int0;
  int i1 = fem::int0;
  int i2m = fem::int0;
  int i2 = fem::int0;
  float q2old = fem::float0;
  float qold = fem::float0;
  float qmov = fem::float0;
  float rbin = fem::float0;
  float rinp = fem::float0;
  float q2new = fem::float0;
  float hc1 = fem::float0;
  float hc2 = fem::float0;
  float ha = fem::float0;
  float pd = fem::float0;
  int im = fem::int0;
  float pes = fem::float0;
  float pqs = fem::float0;
  float fac = fem::float0;
  //C
  //C...Purpose: to modify event so as to approximately take into account
  //C...Bose-Einstein effects according to a simple phenomenological
  //C...parametrization.
  //C
  //C...Boost event to overall CM frame. Calculate CM energy.
  if ((mstj(51) != 1 && mstj(51) != 2) || n - nsav <= 1) {
    return;
  }
  FEM_DO_SAFE(j, 1, 4) {
    dps(j) = 0.e0;
  }
  FEM_DO_SAFE(i, 1, n) {
    if (k(i, 1) <= 0 || k(i, 1) > 10) {
      goto statement_120;
    }
    FEM_DO_SAFE(j, 1, 4) {
      dps(j) += fem::dble(p(i, j));
    }
    statement_120:;
  }
  ludbrb(0, 0, 0.f, 0.f, -dps(1) / dps(4), -dps(2) / dps(4), -dps(3) / dps(4));
  pecm = 0.f;
  FEM_DO_SAFE(i, 1, n) {
    if (k(i, 1) >= 1 && k(i, 1) <= 10) {
      pecm += p(i, 4);
    }
  }
  //C
  //C...Reserve copy of particles by species at end of record.
  nbe(0) = n + mstu(3);
  FEM_DO_SAFE(ibe, 1, fem::min(9, mstj(51))) {
    nbe(ibe) = nbe(ibe - 1);
    FEM_DO_SAFE(i, nsav + 1, n) {
      if (k(i, 2) != kfbe(ibe)) {
        goto statement_150;
      }
      if (k(i, 1) <= 0 || k(i, 1) > 10) {
        goto statement_150;
      }
      if (nbe(ibe) >= mstu(4) - mstu(32) - 5) {
        luerrm(cmn, 11, "(LUBOEI:) no more memory left in LUJETS");
        return;
      }
      nbe(ibe)++;
      k(nbe(ibe), 1) = i;
      FEM_DO_SAFE(j, 1, 3) {
        p(nbe(ibe), j) = 0.f;
      }
      statement_150:;
    }
  }
  //C
  //C...Tabulate integral for subsequent momentum shift.
  FEM_DO_SAFE(ibe, 1, fem::min(9, mstj(51))) {
    if (ibe != 1 && ibe != 4 && ibe <= 7) {
      goto statement_180;
    }
    if (ibe == 1 && fem::max(nbe(1) - nbe(0), nbe(2) - nbe(1), nbe(
        3) - nbe(2)) <= 1) {
      goto statement_180;
    }
    if (ibe == 4 && fem::max(nbe(4) - nbe(3), nbe(5) - nbe(4), nbe(6) - nbe(5),
        nbe(7) - nbe(6)) <= 1) {
      goto statement_180;
    }
    if (ibe >= 8 && nbe(ibe) - nbe(ibe - 1) <= 1) {
      goto statement_180;
    }
    if (ibe == 1) {
      pmhq = 2.f * ulmass(cmn, 211);
    }
    if (ibe == 4) {
      pmhq = 2.f * ulmass(cmn, 321);
    }
    if (ibe == 8) {
      pmhq = 2.f * ulmass(cmn, 221);
    }
    if (ibe == 9) {
      pmhq = 2.f * ulmass(cmn, 331);
    }
    qdel = 0.1f * fem::min(pmhq, parj(93));
    if (mstj(51) == 1) {
      nbin = fem::min(100, fem::nint(9.f * parj(93) / qdel));
      beex = fem::exp(0.5f * qdel / parj(93));
      bert = fem::exp(-qdel / parj(93));
    }
    else {
      nbin = fem::min(100, fem::nint(3.f * parj(93) / qdel));
    }
    FEM_DO_SAFE(ibin, 1, nbin) {
      qbin = qdel * (ibin - 0.5f);
      bei(ibin) = qdel * (fem::pow2(qbin) + fem::pow2(qdel) / 12.f) /
        fem::sqrt(fem::pow2(qbin) + fem::pow2(pmhq));
      if (mstj(51) == 1) {
        beex = beex * bert;
        bei(ibin) = bei(ibin) * beex;
      }
      else {
        bei(ibin) = bei(ibin) * fem::exp(-fem::pow2((qbin / parj(93))));
      }
      if (ibin >= 2) {
        bei(ibin) += bei(ibin - 1);
      }
    }
    //C
    //C...Loop through particle pairs and find old relative momentum.
    statement_180:
    FEM_DO_SAFE(i1m, nbe(ibe - 1) + 1, nbe(ibe) - 1) {
      i1 = k(i1m, 1);
      FEM_DO_SAFE(i2m, i1m + 1, nbe(ibe)) {
        i2 = k(i2m, 1);
        q2old = fem::max(0.f, fem::pow2((p(i1, 4) + p(i2, 4))) -
          fem::pow2((p(i1, 1) + p(i2, 1))) - fem::pow2((p(i1, 2) + p(i2,
          2))) - fem::pow2((p(i1, 3) + p(i2, 3))) - fem::pow2((p(i1,
          5) + p(i2, 5))));
        qold = fem::sqrt(q2old);
        //C
        //C...Calculate new relative momentum.
        if (qold < 0.5f * qdel) {
          qmov = qold / 3.f;
        }
        else if (qold < (nbin - 0.1f) * qdel) {
          rbin = qold / qdel;
          ibin = fem::fint(rbin);
          rinp = (fem::pow3(rbin) - fem::pow3(ibin)) / (3 * ibin * (
            ibin + 1) + 1);
          qmov = (bei(ibin) + rinp * (bei(ibin + 1) - bei(ibin))) *
            fem::sqrt(q2old + fem::pow2(pmhq)) / q2old;
        }
        else {
          qmov = bei(nbin) * fem::sqrt(q2old + fem::pow2(pmhq)) / q2old;
        }
        q2new = q2old * fem::pow((qold / (qold + 3.f * parj(92) * qmov)),
          (2.f / 3.f));
        //C
        //C...Calculate and save shift to be performed on three-momenta.
        hc1 = fem::pow2((p(i1, 4) + p(i2, 4))) - (q2old - q2new);
        hc2 = (q2old - q2new) * fem::pow2((p(i1, 4) - p(i2, 4)));
        ha = 0.5f * (1.f - fem::sqrt(hc1 * q2new / (hc1 * q2old - hc2)));
        FEM_DO_SAFE(j, 1, 3) {
          pd = ha * (p(i2, j) - p(i1, j));
          p(i1m, j) += pd;
          p(i2m, j) = p(i2m, j) - pd;
        }
      }
    }
  }
  //C
  //C...Shift momenta and recalculate energies.
  FEM_DO_SAFE(im, nbe(0) + 1, nbe(fem::min(9, mstj(51)))) {
    i = k(im, 1);
    FEM_DO_SAFE(j, 1, 3) {
      p(i, j) += p(im, j);
    }
    p(i, 4) = fem::sqrt(fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) +
      fem::pow2(p(i, 2)) + fem::pow2(p(i, 3)));
  }
  //C
  //C...Rescale all momenta for energy conservation.
  pes = 0.f;
  pqs = 0.f;
  FEM_DO_SAFE(i, 1, n) {
    if (k(i, 1) <= 0 || k(i, 1) > 10) {
      goto statement_240;
    }
    pes += p(i, 4);
    pqs += fem::pow2(p(i, 5)) / p(i, 4);
    statement_240:;
  }
  fac = (pecm - pqs) / (pes - pqs);
  FEM_DO_SAFE(i, 1, n) {
    if (k(i, 1) <= 0 || k(i, 1) > 10) {
      goto statement_260;
    }
    FEM_DO_SAFE(j, 1, 3) {
      p(i, j) = fac * p(i, j);
    }
    p(i, 4) = fem::sqrt(fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) +
      fem::pow2(p(i, 2)) + fem::pow2(p(i, 3)));
    statement_260:;
  }
  //C
  //C...Boost back to correct reference frame.
  ludbrb(0, 0, 0.f, 0.f, dps(1) / dps(4), dps(2) / dps(4), dps(3) / dps(4));
  //C
}

//C
//C*********************************************************************
//C
void
luexec(
  common& cmn)
{
  int& n = cmn.n;
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_cref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  //
  int mcons = fem::int0;
  int nsav = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  arr_2d<2, 6, float> ps(fem::fill0);
  int mbe = fem::int0;
  int ip = fem::int0;
  int kc = fem::int0;
  int ip1 = fem::int0;
  float qmax = fem::float0;
  float pip5 = fem::float0;
  int mfrag = fem::int0;
  float pdev = fem::float0;
  //C
  //C...Purpose: to administrate the fragmentation and decay chain.
  //C
  //C...Initialize and reset.
  mstu(24) = 0;
  if (mstu(12) >= 1) {
    lulist(cmn, 0);
  }
  mstu(31)++;
  mstu(1) = 0;
  mstu(2) = 0;
  mstu(3) = 0;
  mcons = 1;
  //C
  //C...Sum up momentum, energy and charge for starting entries.
  nsav = n;
  FEM_DO_SAFE(i, 1, 2) {
    FEM_DO_SAFE(j, 1, 6) {
      ps(i, j) = 0.f;
    }
  }
  FEM_DO_SAFE(i, 1, n) {
    if (k(i, 1) <= 0 || k(i, 1) > 10) {
      goto statement_120;
    }
    FEM_DO_SAFE(j, 1, 4) {
      ps(1, j) += p(i, j);
    }
    ps(1, 6) += luchge(cmn, k(i, 2));
    statement_120:;
  }
  paru(21) = ps(1, 4);
  //C
  //C...Prepare system for subsequent fragmentation/decay.
  luprep(cmn, 0);
  //C
  //C...Loop through jet fragmentation and particle decays.
  mbe = 0;
  statement_130:
  mbe++;
  ip = 0;
  statement_140:
  ip++;
  kc = 0;
  if (k(ip, 1) > 0 && k(ip, 1) <= 10) {
    kc = lucomp(cmn, k(ip, 2));
  }
  if (kc == 0) {
    //C
    //C...Particle decay if unstable and allowed. Save long-lived particle
    //C...decays until second pass after Bose-Einstein effects.
  }
  else if (kchg(kc, 2) == 0) {
    //Clin-4/2008 break up compound IF statements:
    //C        IF(MSTJ(21).GE.1.AND.MDCY(KC,1).GE.1.AND.(MSTJ(51).LE.0.OR.MBE.
    //C     &  EQ.2.OR.PMAS(KC,2).GE.PARJ(91).OR.IABS(K(IP,2)).EQ.311))
    //C     &  CALL LUDECY(IP)
    if (mstj(21) >= 1 && mdcy(kc, 1) >= 1) {
      if (mstj(51) <= 0 || mbe == 2 || pmas(kc, 2) >= parj(91) ||
          fem::iabs(k(ip, 2)) == 311) {
        ludecy(cmn, ip);
      }
    }
    //C
    //C...Decay products may develop a shower.
    if (mstj(92) > 0) {
      ip1 = mstj(92);
      qmax = fem::sqrt(fem::max(0.f, fem::pow2((p(ip1, 4) + p(ip1 + 1,
        4))) - fem::pow2((p(ip1, 1) + p(ip1 + 1, 1))) - fem::pow2((p(ip1,
        2) + p(ip1 + 1, 2))) - fem::pow2((p(ip1, 3) + p(ip1 + 1,
        3)))));
      lushow(cmn, ip1, ip1 + 1, qmax);
      luprep(cmn, ip1);
      mstj(92) = 0;
    }
    else if (mstj(92) < 0) {
      ip1 = -mstj(92);
      //Clin-8/19/02 avoid actual argument in common blocks of LUSHOW:
      //C          CALL LUSHOW(IP1,-3,P(IP,5))
      pip5 = p(ip, 5);
      lushow(cmn, ip1, -3, pip5);
      luprep(cmn, ip1);
      mstj(92) = 0;
    }
    //C
    //C...Jet fragmentation: string or independent fragmentation.
  }
  else if (k(ip, 1) == 1 || k(ip, 1) == 2) {
    mfrag = mstj(1);
    if (mfrag >= 1 && k(ip, 1) == 1) {
      mfrag = 2;
    }
    if (mstj(21) >= 2 && k(ip, 1) == 2 && n > ip) {
      if (k(ip + 1, 1) == 1 && k(ip + 1, 3) == k(ip, 3) && k(ip,
          3) > 0 && k(ip, 3) < ip) {
        if (kchg(lucomp(cmn, k(k(ip, 3), 2)), 2) == 0) {
          mfrag = fem::min(1, mfrag);
        }
      }
    }
    if (mfrag == 1) {
      lustrf(cmn, ip);
    }
    if (mfrag == 2) {
      luindf(cmn, ip);
    }
    if (mfrag == 2 && k(ip, 1) == 1) {
      mcons = 0;
    }
    if (mfrag == 2 && (mstj(3) <= 0 || fem::mod(mstj(3), 5) == 0)) {
      mcons = 0;
    }
  }
  //C
  //C...Loop back if enough space left in LUJETS and no error abort.
  if (mstu(24) != 0 && mstu(21) >= 2) {
  }
  else if (ip < n && n < mstu(4) - 20 - mstu(32)) {
    goto statement_140;
  }
  else if (ip < n) {
    luerrm(cmn, 11, "(LUEXEC:) no more memory left in LUJETS");
  }
  //C
  //C...Include simple Bose-Einstein effect parametrization if desired.
  if (mbe == 1 && mstj(51) >= 1) {
    luboei(cmn, nsav);
    goto statement_130;
  }
  //C
  //C...Check that momentum, energy and charge were conserved.
  FEM_DO_SAFE(i, 1, n) {
    if (k(i, 1) <= 0 || k(i, 1) > 10) {
      goto statement_160;
    }
    FEM_DO_SAFE(j, 1, 4) {
      ps(2, j) += p(i, j);
    }
    ps(2, 6) += luchge(cmn, k(i, 2));
    statement_160:;
  }
  pdev = (fem::abs(ps(2, 1) - ps(1, 1)) + fem::abs(ps(2, 2) - ps(1,
    2)) + fem::abs(ps(2, 3) - ps(1, 3)) + fem::abs(ps(2, 4) - ps(1,
    4))) / (1.f + fem::abs(ps(2, 4)) + fem::abs(ps(1, 4)));
  if (mcons == 1 && pdev > paru(11)) {
    luerrm(cmn, 15, "(LUEXEC:) four-momentum was not conserved");
  }
  //C      IF(MCONS.EQ.1.AND.PDEV.GT.PARU(11)) then
  //C         CALL LUERRM(15,
  //C     &'(LUEXEC:) four-momentum was not conserved')
  //C         write(6,*) 'PS1,2=',PS(1,1),PS(1,2),PS(1,3),PS(1,4),
  //C     1        '*',PS(2,1),PS(2,2),PS(2,3),PS(2,4)
  //C      endif
  //C
  if (mcons == 1 && fem::abs(ps(2, 6) - ps(1, 6)) > 0.1f) {
    luerrm(cmn, 15, "(LUEXEC:) charge was not conserved");
  }
  //C
}

//C.................... hipyset1.35.f
//C
//C     Modified for HIJING program
//C
//C    modification July 22, 1997  In pyremnn put an upper limit
//C     on the total pt kick the parton can accumulate via multiple
//C     scattering. Set the upper limit to be the sqrt(s)/2,
//C     this is fix cronin bug for Pb+Pb events at SPS energy.
//C
//C Last modification Oct. 1993 to comply with non-vax
//C machines' compiler
//C
//C*********************************************************************
//C
void
lu2ent(
  common& cmn,
  int const& ip,
  int const& kf1,
  int const& kf2,
  float const& pecm)
{
  // COMMON lujets
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  // COMMON ludat1
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  // COMMON ludat2
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  //C
  //C...Purpose: to store two partons/particles in their CM frame,
  //C...with the first along the +z axis.
  //C
  //C...Standard checks.
  mstu(28) = 0;
  if (mstu(12) >= 1) {
    lulist(cmn, 0);
  }
  int ipa = fem::max(1, fem::iabs(ip));
  if (ipa > mstu(4) - 1) {
    luerrm(cmn, 21, "(LU2ENT:) writing outside LUJETS memory");
  }
  int kc1 = lucomp(cmn, kf1);
  int kc2 = lucomp(cmn, kf2);
  if (kc1 == 0 || kc2 == 0) {
    luerrm(cmn, 12, "(LU2ENT:) unknown flavour code");
  }
  //C
  //C...Find masses. Reset K, P and V vectors.
  float pm1 = 0.f;
  if (mstu(10) == 1) {
    pm1 = p(ipa, 5);
  }
  if (mstu(10) >= 2) {
    pm1 = ulmass(cmn, kf1);
  }
  float pm2 = 0.f;
  if (mstu(10) == 1) {
    pm2 = p(ipa + 1, 5);
  }
  if (mstu(10) >= 2) {
    pm2 = ulmass(cmn, kf2);
  }
  int i = fem::int0;
  int j = fem::int0;
  FEM_DO_SAFE(i, ipa, ipa + 1) {
    FEM_DO_SAFE(j, 1, 5) {
      k(i, j) = 0;
      p(i, j) = 0.f;
      v(i, j) = 0.f;
    }
  }
  //C
  //C...Check flavours.
  int kq1 = kchg(kc1, 2) * fem::isign(1, kf1);
  int kq2 = kchg(kc2, 2) * fem::isign(1, kf2);
  if (kq1 + kq2 != 0 && kq1 + kq2 != 4) {
    luerrm(cmn, 2, "(LU2ENT:) unphysical flavour combination");
  }
  k(ipa, 2) = kf1;
  k(ipa + 1, 2) = kf2;
  //C
  //C...Store partons/particles in K vectors for normal case.
  if (ip >= 0) {
    k(ipa, 1) = 1;
    if (kq1 != 0 && kq2 != 0) {
      k(ipa, 1) = 2;
    }
    k(ipa + 1, 1) = 1;
    //C
    //C...Store partons in K vectors for parton shower evolution.
  }
  else {
    if (kq1 == 0 || kq2 == 0) {
      luerrm(cmn, 2,
        "(LU2ENT:) requested flavours can not develop parton shower");
    }
    k(ipa, 1) = 3;
    k(ipa + 1, 1) = 3;
    k(ipa, 4) = mstu(5) * (ipa + 1);
    k(ipa, 5) = k(ipa, 4);
    k(ipa + 1, 4) = mstu(5) * ipa;
    k(ipa + 1, 5) = k(ipa + 1, 4);
  }
  //C
  //C...Check kinematics and store partons/particles in P vectors.
  if (pecm <= pm1 + pm2) {
    luerrm(cmn, 13, "(LU2ENT:) energy smaller than sum of masses");
  }
  float pa = fem::sqrt(fem::max(0.f, fem::pow2((fem::pow2(pecm) -
    fem::pow2(pm1) - fem::pow2(pm2))) - fem::pow2((2.f * pm1 * pm2)))) /
    (2.f * pecm);
  p(ipa, 3) = pa;
  p(ipa, 4) = fem::sqrt(fem::pow2(pm1) + fem::pow2(pa));
  p(ipa, 5) = pm1;
  p(ipa + 1, 3) = -pa;
  p(ipa + 1, 4) = fem::sqrt(fem::pow2(pm2) + fem::pow2(pa));
  p(ipa + 1, 5) = pm2;
  //C
  //C...Set N. Optionally fragment/decay.
  cmn.n = ipa + 1;
  if (ip == 0) {
    luexec(cmn);
  }
  //C
}

struct lugive_save
{
  arr<fem::str<26> > chalp;
  arr<fem::str<4> > chvar;

  lugive_save() :
    chalp(dimension(2), fem::fill0),
    chvar(dimension(17), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
lugive(
  common& cmn,
  str_cref chin)
{
  FEM_CMN_SVE(lugive);
  common_read read(cmn);
  common_write write(cmn);
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_ref<float> parj(cmn.parj, dimension(200));
  arr_ref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_ref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_ref<float> parf(cmn.parf, dimension(2000));
  arr_ref<float, 2> vckm(cmn.vckm, dimension(4, 4));
  arr_ref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_ref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_ref<float> brat(cmn.brat, dimension(2000));
  arr_ref<int, 2> kfdp(cmn.kfdp, dimension(2000, 5));
  str_arr_ref<1> chaf(cmn.chaf, dimension(500));
  //
  str_arr_ref<1> chalp(sve.chalp, dimension(2));
  str_arr_ref<1> chvar(sve.chvar, dimension(17));
  if (is_called_first_time) {
    {
      static const char* values[] = {
        "N", "K", "P", "V", "MSTU", "PARU", "MSTJ", "PARJ", "KCHG",
          "PMAS", "PARF", "VCKM", "MDCY", "MDME", "BRAT", "KFDP",
          "CHAF"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chvar;
    }
    {
      static const char* values[] = {
        "abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chalp;
    }
  }
  fem::str<104> chbit = fem::char0;
  int lbit = fem::int0;
  int ltot = fem::int0;
  int lcom = fem::int0;
  fem::str<104> chfix = fem::char0;
  int llow = fem::int0;
  int lhig = fem::int0;
  int lnam = fem::int0;
  fem::str<4> chnam = fem::char0;
  int lalp = fem::int0;
  int ivar = fem::int0;
  int iv = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  int lind = fem::int0;
  fem::str<8> chind = fem::char0;
  int i1 = fem::int0;
  int ierr = fem::int0;
  int iold = fem::int0;
  float rold = fem::float0;
  fem::str<8> chold = fem::char0;
  fem::str<10> chini = fem::char0;
  int inew = fem::int0;
  fem::str<16> chinr = fem::char0;
  float rnew = fem::float0;
  fem::str<8> chnew = fem::char0;
  static const char* format_1000 = "(5x,a60)";
  //C
  //C...Purpose: to set values of commonblock variables.
  //C
  //C...Length of character variable. Subdivide it into instructions.
  if (mstu(12) >= 1) {
    lulist(cmn, 0);
  }
  chbit = chin + str_cref(" ");
  lbit = 101;
  statement_100:
  lbit = lbit - 1;
  if (chbit(lbit, lbit) == " ") {
    goto statement_100;
  }
  ltot = 0;
  FEM_DO_SAFE(lcom, 1, lbit) {
    if (chbit(lcom, lcom) == " ") {
      goto statement_110;
    }
    ltot++;
    chfix(ltot, ltot) = chbit(lcom, lcom);
    statement_110:;
  }
  llow = 0;
  statement_120:
  lhig = llow + 1;
  statement_130:
  lhig++;
  if (lhig <= ltot && chfix(lhig, lhig) != ";") {
    goto statement_130;
  }
  lbit = lhig - llow - 1;
  chbit(1, lbit) = chfix(llow + 1, lhig - 1);
  //C
  //C...Identify commonblock variable.
  lnam = 1;
  statement_140:
  lnam++;
  if (chbit(lnam, lnam) != "(" && chbit(lnam, lnam) != "=" && lnam <= 4) {
    goto statement_140;
  }
  chnam = chbit(1, lnam - 1) + str_cref(" ");
  FEM_DO_SAFE(lcom, 1, lnam - 1) {
    FEM_DO_SAFE(lalp, 1, 26) {
      if (chnam(lcom, lcom) == chalp(1)(lalp, lalp)) {
        chnam(lcom, lcom) = chalp(2)(lalp, lalp);
      }
    }
  }
  ivar = 0;
  FEM_DO_SAFE(iv, 1, 17) {
    if (chnam == chvar(iv)) {
      ivar = iv;
    }
  }
  if (ivar == 0) {
    luerrm(cmn, 18, "(LUGIVE:) do not recognize variable " + chnam);
    llow = lhig;
    if (llow < ltot) {
      goto statement_120;
    }
    return;
  }
  //C
  //C...Identify any indices.
  i = 0;
  j = 0;
  if (chbit(lnam, lnam) == "(") {
    lind = lnam;
    statement_170:
    lind++;
    if (chbit(lind, lind) != ")" && chbit(lind, lind) != ",") {
      goto statement_170;
    }
    chind = " ";
    if ((chbit(lnam + 1, lnam + 1) == "C" || chbit(lnam + 1, lnam +
        1) == "c") && (ivar == 9 || ivar == 10 || ivar == 13 ||
        ivar == 17)) {
      chind(lnam - lind + 11, 8) = chbit(lnam + 2, lind - 1);
      read(chind, "(i8)"), i1;
      i = lucomp(cmn, i1);
    }
    else {
      chind(lnam - lind + 10, 8) = chbit(lnam + 1, lind - 1);
      read(chind, "(i8)"), i;
    }
    lnam = lind;
    if (chbit(lnam, lnam) == ")") {
      lnam++;
    }
  }
  if (chbit(lnam, lnam) == ",") {
    lind = lnam;
    statement_180:
    lind++;
    if (chbit(lind, lind) != ")" && chbit(lind, lind) != ",") {
      goto statement_180;
    }
    chind = " ";
    chind(lnam - lind + 10, 8) = chbit(lnam + 1, lind - 1);
    read(chind, "(i8)"), j;
    lnam = lind + 1;
  }
  //C
  //C...Check that indices allowed and save old value.
  ierr = 1;
  if (chbit(lnam, lnam) != "=") {
    goto statement_190;
  }
  if (ivar == 1) {
    if (i != 0 || j != 0) {
      goto statement_190;
    }
    iold = n;
  }
  else if (ivar == 2) {
    if (i < 1 || i > mstu(4) || j < 1 || j > 5) {
      goto statement_190;
    }
    iold = k(i, j);
  }
  else if (ivar == 3) {
    if (i < 1 || i > mstu(4) || j < 1 || j > 5) {
      goto statement_190;
    }
    rold = p(i, j);
  }
  else if (ivar == 4) {
    if (i < 1 || i > mstu(4) || j < 1 || j > 5) {
      goto statement_190;
    }
    rold = v(i, j);
  }
  else if (ivar == 5) {
    if (i < 1 || i > 200 || j != 0) {
      goto statement_190;
    }
    iold = mstu(i);
  }
  else if (ivar == 6) {
    if (i < 1 || i > 200 || j != 0) {
      goto statement_190;
    }
    rold = paru(i);
  }
  else if (ivar == 7) {
    if (i < 1 || i > 200 || j != 0) {
      goto statement_190;
    }
    iold = mstj(i);
  }
  else if (ivar == 8) {
    if (i < 1 || i > 200 || j != 0) {
      goto statement_190;
    }
    rold = parj(i);
  }
  else if (ivar == 9) {
    if (i < 1 || i > mstu(6) || j < 1 || j > 3) {
      goto statement_190;
    }
    iold = kchg(i, j);
  }
  else if (ivar == 10) {
    if (i < 1 || i > mstu(6) || j < 1 || j > 4) {
      goto statement_190;
    }
    rold = pmas(i, j);
  }
  else if (ivar == 11) {
    if (i < 1 || i > 2000 || j != 0) {
      goto statement_190;
    }
    rold = parf(i);
  }
  else if (ivar == 12) {
    if (i < 1 || i > 4 || j < 1 || j > 4) {
      goto statement_190;
    }
    rold = vckm(i, j);
  }
  else if (ivar == 13) {
    if (i < 1 || i > mstu(6) || j < 1 || j > 3) {
      goto statement_190;
    }
    iold = mdcy(i, j);
  }
  else if (ivar == 14) {
    if (i < 1 || i > mstu(7) || j < 1 || j > 2) {
      goto statement_190;
    }
    iold = mdme(i, j);
  }
  else if (ivar == 15) {
    if (i < 1 || i > mstu(7) || j != 0) {
      goto statement_190;
    }
    rold = brat(i);
  }
  else if (ivar == 16) {
    if (i < 1 || i > mstu(7) || j < 1 || j > 5) {
      goto statement_190;
    }
    iold = kfdp(i, j);
  }
  else if (ivar == 17) {
    if (i < 1 || i > mstu(6) || j != 0) {
      goto statement_190;
    }
    chold = chaf(i);
  }
  ierr = 0;
  statement_190:
  if (ierr == 1) {
    luerrm(cmn, 18, "(LUGIVE:) unallowed indices for " + chbit(1, lnam - 1));
    llow = lhig;
    if (llow < ltot) {
      goto statement_120;
    }
    return;
  }
  //C
  //C...Print current value of variable. Loop back.
  if (lnam >= lbit) {
    chbit(lnam, 14) = " ";
    chbit(15, 60) = " has the value                                ";
    if (ivar == 1 || ivar == 2 || ivar == 5 || ivar == 7 ||
        ivar == 9 || ivar == 13 || ivar == 14 || ivar == 16) {
      write(chbit(51, 60), "(i10)"), iold;
    }
    else if (ivar != 17) {
      write(chbit(47, 60), "(f14.5)"), rold;
    }
    else {
      chbit(53, 60) = chold;
    }
    if (mstu(13) >= 1) {
      write(mstu(11), format_1000), chbit(1, 60);
    }
    llow = lhig;
    if (llow < ltot) {
      goto statement_120;
    }
    return;
  }
  //C
  //C...Read in new variable value.
  if (ivar == 1 || ivar == 2 || ivar == 5 || ivar == 7 ||
      ivar == 9 || ivar == 13 || ivar == 14 || ivar == 16) {
    chini = " ";
    chini(lnam - lbit + 11, 10) = chbit(lnam + 1, lbit);
    read(chini, "(i10)"), inew;
  }
  else if (ivar != 17) {
    chinr = " ";
    chinr(lnam - lbit + 17, 16) = chbit(lnam + 1, lbit);
    read(chinr, "(f16.2)"), rnew;
  }
  else {
    chnew = chbit(lnam + 1, lbit) + str_cref(" ");
  }
  //C
  //C...Store new variable value.
  if (ivar == 1) {
    n = inew;
  }
  else if (ivar == 2) {
    k(i, j) = inew;
  }
  else if (ivar == 3) {
    p(i, j) = rnew;
  }
  else if (ivar == 4) {
    v(i, j) = rnew;
  }
  else if (ivar == 5) {
    mstu(i) = inew;
  }
  else if (ivar == 6) {
    paru(i) = rnew;
  }
  else if (ivar == 7) {
    mstj(i) = inew;
  }
  else if (ivar == 8) {
    parj(i) = rnew;
  }
  else if (ivar == 9) {
    kchg(i, j) = inew;
  }
  else if (ivar == 10) {
    pmas(i, j) = rnew;
  }
  else if (ivar == 11) {
    parf(i) = rnew;
  }
  else if (ivar == 12) {
    vckm(i, j) = rnew;
  }
  else if (ivar == 13) {
    mdcy(i, j) = inew;
  }
  else if (ivar == 14) {
    mdme(i, j) = inew;
  }
  else if (ivar == 15) {
    brat(i) = rnew;
  }
  else if (ivar == 16) {
    kfdp(i, j) = inew;
  }
  else if (ivar == 17) {
    chaf(i) = chnew;
  }
  //C
  //C...Write old and new value. Loop back.
  chbit(lnam, 14) = " ";
  chbit(15, 60) = " changed from                to               ";
  if (ivar == 1 || ivar == 2 || ivar == 5 || ivar == 7 ||
      ivar == 9 || ivar == 13 || ivar == 14 || ivar == 16) {
    write(chbit(33, 42), "(i10)"), iold;
    write(chbit(51, 60), "(i10)"), inew;
  }
  else if (ivar != 17) {
    write(chbit(29, 42), "(f14.5)"), rold;
    write(chbit(47, 60), "(f14.5)"), rnew;
  }
  else {
    chbit(35, 42) = chold;
    chbit(53, 60) = chnew;
  }
  if (mstu(13) >= 1) {
    write(mstu(11), format_1000), chbit(1, 60);
  }
  llow = lhig;
  if (llow < ltot) {
    goto statement_120;
  }
  //C
  //C...Format statement for output on unit MSTU(11) (by default 6).
  //C
}

//C
//C*********************************************************************
//C
float
ulalps(
  common& cmn,
  float const& q2)
{
  float return_value = fem::float0;
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<float> paru(cmn.paru, dimension(200));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  //
  float q2eff = fem::float0;
  int nf = fem::int0;
  float alam2 = fem::float0;
  float q2thr = fem::float0;
  float b0 = fem::float0;
  float algq = fem::float0;
  float b1 = fem::float0;
  //C
  //C...Purpose: to give the value of alpha_strong.
  //C
  //C...Constant alpha_strong trivial.
  if (mstu(111) <= 0) {
    return_value = paru(111);
    mstu(118) = mstu(112);
    paru(117) = 0.f;
    paru(118) = paru(111);
    return return_value;
  }
  //C
  //C...Find effective Q2, number of flavours and Lambda.
  q2eff = q2;
  if (mstu(115) >= 2) {
    q2eff = fem::max(q2, paru(114));
  }
  nf = mstu(112);
  alam2 = fem::pow2(paru(112));
  statement_100:
  if (nf > fem::max(2, mstu(113))) {
    q2thr = paru(113) * fem::pow2(pmas(nf, 1));
    if (q2eff < q2thr) {
      nf = nf - 1;
      alam2 = alam2 * fem::pow((q2thr / alam2), (2.f / (33.f - 2.f * nf)));
      goto statement_100;
    }
  }
  statement_110:
  if (nf < fem::min(8, mstu(114))) {
    q2thr = paru(113) * fem::pow2(pmas(nf + 1, 1));
    if (q2eff > q2thr) {
      nf++;
      alam2 = alam2 * fem::pow((alam2 / q2thr), (2.f / (33.f - 2.f * nf)));
      goto statement_110;
    }
  }
  if (mstu(115) == 1) {
    q2eff += alam2;
  }
  paru(117) = fem::sqrt(alam2);
  //C
  //C...Evaluate first or second order alpha_strong.
  b0 = (33.f - 2.f * nf) / 6.f;
  algq = fem::log(q2eff / alam2);
  if (mstu(111) == 1) {
    return_value = paru(2) / (b0 * algq);
  }
  else {
    b1 = (153.f - 19.f * nf) / 6.f;
    return_value = paru(2) / (b0 * algq) * (1.f - b1 * fem::log(
      algq) / (fem::pow2(b0) * algq));
  }
  mstu(118) = nf;
  paru(118) = return_value;
  //C
  return return_value;
}

//C
//C*********************************************************************
//C
void
lurobo(
  common& cmn,
  float const& the,
  float const& phi,
  float const& bex,
  float const& bey,
  float const& bez)
{
  int& n = cmn.n;
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  //
  int imin = fem::int0;
  int imax = fem::int0;
  double dbx = fem::double0;
  double dby = fem::double0;
  double dbz = fem::double0;
  int imi = fem::int0;
  int ima = fem::int0;
  double dbex = fem::double0;
  double dbey = fem::double0;
  double dbez = fem::double0;
  arr_2d<3, 3, float> rot(fem::fill0);
  int i = fem::int0;
  int j = fem::int0;
  arr_1d<3, float> pr(fem::fill0);
  arr_1d<3, float> vr(fem::fill0);
  double db = fem::double0;
  double dga = fem::double0;
  arr_1d<4, double> dp(fem::fill0);
  arr_1d<4, double> dv(fem::fill0);
  double dbp = fem::double0;
  double dgabp = fem::double0;
  double dbv = fem::double0;
  double dgabv = fem::double0;
  //C
  //C...Purpose: to perform rotations and boosts.
  //C
  //C...Find range of rotation/boost. Convert boost to double precision.
  imin = 1;
  if (mstu(1) > 0) {
    imin = mstu(1);
  }
  imax = n;
  if (mstu(2) > 0) {
    imax = mstu(2);
  }
  dbx = fem::dble(bex);
  dby = fem::dble(bey);
  dbz = fem::dble(bez);
  goto statement_100;
  //C
  //C...Entry for specific range and double precision boost.
  // UNHANDLED: ENTRY ludbrb(imi,ima,the,phi,dbex,dbey,dbez)
  imin = imi;
  if (imin <= 0) {
    imin = 1;
  }
  imax = ima;
  if (imax <= 0) {
    imax = n;
  }
  dbx = dbex;
  dby = dbey;
  dbz = dbez;
  //C
  //C...Check range of rotation/boost.
  statement_100:
  if (imin > mstu(4) || imax > mstu(4)) {
    luerrm(cmn, 11, "(LUROBO:) range outside LUJETS memory");
    return;
  }
  //C
  //C...Rotate, typically from z axis to direction (theta,phi).
  //Clin-5/2012:
  //C      IF(THE**2+PHI**2.GT.1E-20) THEN
  if ((fem::pow2(the) + fem::pow2(phi)) > 1e-20f) {
    rot(1, 1) = fem::cos(the) * fem::cos(phi);
    rot(1, 2) = -fem::sin(phi);
    rot(1, 3) = fem::sin(the) * fem::cos(phi);
    rot(2, 1) = fem::cos(the) * fem::sin(phi);
    rot(2, 2) = fem::cos(phi);
    rot(2, 3) = fem::sin(the) * fem::sin(phi);
    rot(3, 1) = -fem::sin(the);
    rot(3, 2) = 0.f;
    rot(3, 3) = fem::cos(the);
    FEM_DO_SAFE(i, imin, imax) {
      if (k(i, 1) <= 0) {
        goto statement_130;
      }
      FEM_DO_SAFE(j, 1, 3) {
        pr(j) = p(i, j);
        vr(j) = v(i, j);
      }
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) = rot(j, 1) * pr(1) + rot(j, 2) * pr(2) + rot(j, 3) * pr(3);
        v(i, j) = rot(j, 1) * vr(1) + rot(j, 2) * vr(2) + rot(j, 3) * vr(3);
      }
      statement_130:;
    }
  }
  //C
  //C...Boost, typically from rest to momentum/energy=beta.
  //Clin-5/2012:
  //C      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
  if ((fem::pow2(dbx) + fem::pow2(dby) + fem::pow2(dbz)) > 1e-20) {
    db = fem::sqrt(fem::pow2(dbx) + fem::pow2(dby) + fem::pow2(dbz));
    if (db > 0.99999999e0) {
      //C...Rescale boost vector if too close to unity.
      luerrm(cmn, 3, "(LUROBO:) boost vector too large");
      dbx = dbx * (0.99999999e0 / db);
      dby = dby * (0.99999999e0 / db);
      dbz = dbz * (0.99999999e0 / db);
      db = 0.99999999e0;
    }
    dga = 1e0 / fem::sqrt(1e0 - fem::pow2(db));
    FEM_DO_SAFE(i, imin, imax) {
      if (k(i, 1) <= 0) {
        goto statement_150;
      }
      FEM_DO_SAFE(j, 1, 4) {
        dp(j) = fem::dble(p(i, j));
        dv(j) = fem::dble(v(i, j));
      }
      dbp = dbx * dp(1) + dby * dp(2) + dbz * dp(3);
      dgabp = dga * (dga * dbp / (1e0 + dga) + dp(4));
      p(i, 1) = fem::sngl(dp(1) + dgabp * dbx);
      p(i, 2) = fem::sngl(dp(2) + dgabp * dby);
      p(i, 3) = fem::sngl(dp(3) + dgabp * dbz);
      p(i, 4) = fem::sngl(dga * (dp(4) + dbp));
      dbv = dbx * dv(1) + dby * dv(2) + dbz * dv(3);
      dgabv = dga * (dga * dbv / (1e0 + dga) + dv(4));
      v(i, 1) = fem::sngl(dv(1) + dgabv * dbx);
      v(i, 2) = fem::sngl(dv(2) + dgabv * dby);
      v(i, 3) = fem::sngl(dv(3) + dgabv * dbz);
      v(i, 4) = fem::sngl(dga * (dv(4) + dbv));
      statement_150:;
    }
  }
  //C
}

//C
//C*********************************************************************
//C THIS SUBROUTINE IS ONLY FOR THE USE OF HIJING TO ROTATE OR BOOST
//C        THE FOUR MOMENTUM ONLY
//C*********************************************************************
//C
void
hirobo(
  common& cmn,
  float const& the,
  float const& phi,
  float const& bex,
  float const& bey,
  float const& bez)
{
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  //
  int imin = fem::int0;
  int imax = fem::int0;
  double dbx = fem::double0;
  double dby = fem::double0;
  double dbz = fem::double0;
  arr_2d<3, 3, float> rot(fem::fill0);
  int i = fem::int0;
  int j = fem::int0;
  arr_1d<3, float> pr(fem::fill0);
  double db = fem::double0;
  double dga = fem::double0;
  arr_1d<4, double> dp(fem::fill0);
  double dbp = fem::double0;
  double dgabp = fem::double0;
  //C
  //C...Purpose: to perform rotations and boosts.
  //C
  //C...Find range of rotation/boost. Convert boost to double precision.
  imin = 1;
  if (mstu(1) > 0) {
    imin = mstu(1);
  }
  imax = cmn.n;
  if (mstu(2) > 0) {
    imax = mstu(2);
  }
  dbx = fem::dble(bex);
  dby = fem::dble(bey);
  dbz = fem::dble(bez);
  //C
  //C...Check range of rotation/boost.
  if (imin > mstu(4) || imax > mstu(4)) {
    luerrm(cmn, 11, "(LUROBO:) range outside LUJETS memory");
    return;
  }
  //C
  //C...Rotate, typically from z axis to direction (theta,phi).
  //Clin-5/2012:
  //C      IF(THE**2+PHI**2.GT.1E-20) THEN
  if ((fem::pow2(the) + fem::pow2(phi)) > 1e-20f) {
    rot(1, 1) = fem::cos(the) * fem::cos(phi);
    rot(1, 2) = -fem::sin(phi);
    rot(1, 3) = fem::sin(the) * fem::cos(phi);
    rot(2, 1) = fem::cos(the) * fem::sin(phi);
    rot(2, 2) = fem::cos(phi);
    rot(2, 3) = fem::sin(the) * fem::sin(phi);
    rot(3, 1) = -fem::sin(the);
    rot(3, 2) = 0.f;
    rot(3, 3) = fem::cos(the);
    FEM_DO_SAFE(i, imin, imax) {
      if (k(i, 1) <= 0) {
        goto statement_130;
      }
      FEM_DO_SAFE(j, 1, 3) {
        pr(j) = p(i, j);
      }
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) = rot(j, 1) * pr(1) + rot(j, 2) * pr(2) + rot(j, 3) * pr(3);
      }
      statement_130:;
    }
  }
  //C
  //C...Boost, typically from rest to momentum/energy=beta.
  //Clin-5/2012:
  //C      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
  if ((fem::pow2(dbx) + fem::pow2(dby) + fem::pow2(dbz)) > 1e-20) {
    db = fem::sqrt(fem::pow2(dbx) + fem::pow2(dby) + fem::pow2(dbz));
    if (db > 0.99999999e0) {
      //C...Rescale boost vector if too close to unity.
      luerrm(cmn, 3, "(LUROBO:) boost vector too large");
      dbx = dbx * (0.99999999e0 / db);
      dby = dby * (0.99999999e0 / db);
      dbz = dbz * (0.99999999e0 / db);
      db = 0.99999999e0;
    }
    dga = 1e0 / fem::sqrt(1e0 - fem::pow2(db));
    FEM_DO_SAFE(i, imin, imax) {
      if (k(i, 1) <= 0) {
        goto statement_150;
      }
      FEM_DO_SAFE(j, 1, 4) {
        dp(j) = fem::dble(p(i, j));
      }
      dbp = dbx * dp(1) + dby * dp(2) + dbz * dp(3);
      dgabp = dga * (dga * dbp / (1e0 + dga) + dp(4));
      p(i, 1) = fem::sngl(dp(1) + dgabp * dbx);
      p(i, 2) = fem::sngl(dp(2) + dgabp * dby);
      p(i, 3) = fem::sngl(dp(3) + dgabp * dbz);
      p(i, 4) = fem::sngl(dga * (dp(4) + dbp));
      statement_150:;
    }
  }
  //C
}

//C
//C*********************************************************************
//C
void
luedit(
  common& cmn,
  int const& medit)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  //
  int imax = fem::int0;
  int i1 = fem::int0;
  int i = fem::int0;
  int kc = fem::int0;
  int j = fem::int0;
  int id = fem::int0;
  int im = fem::int0;
  int kcm = fem::int0;
  int kcd = fem::int0;
  int kh = fem::int0;
  int is = fem::int0;
  arr_1d<2, int> ns(fem::fill0);
  arr_1d<2, float> pts(fem::fill0);
  arr_1d<2, float> pls(fem::fill0);
  //C
  //C...Purpose: to perform global manipulations on the event record,
  //C...in particular to exclude unstable or undetectable partons/particles.
  //C
  //C...Remove unwanted partons/particles.
  if ((medit >= 0 && medit <= 3) || medit == 5) {
    imax = n;
    if (mstu(2) > 0) {
      imax = mstu(2);
    }
    i1 = fem::max(1, mstu(1)) - 1;
    FEM_DO_SAFE(i, fem::max(1, mstu(1)), imax) {
      if (k(i, 1) == 0 || k(i, 1) > 20) {
        goto statement_110;
      }
      if (medit == 1) {
        if (k(i, 1) > 10) {
          goto statement_110;
        }
      }
      else if (medit == 2) {
        if (k(i, 1) > 10) {
          goto statement_110;
        }
        kc = lucomp(cmn, k(i, 2));
        if (kc == 0 || kc == 12 || kc == 14 || kc == 16 || kc == 18) {
          goto statement_110;
        }
      }
      else if (medit == 3) {
        if (k(i, 1) > 10) {
          goto statement_110;
        }
        kc = lucomp(cmn, k(i, 2));
        if (kc == 0) {
          goto statement_110;
        }
        if (kchg(kc, 2) == 0 && luchge(cmn, k(i, 2)) == 0) {
          goto statement_110;
        }
      }
      else if (medit == 5) {
        if (k(i, 1) == 13 || k(i, 1) == 14) {
          goto statement_110;
        }
        kc = lucomp(cmn, k(i, 2));
        if (kc == 0) {
          goto statement_110;
        }
        if (k(i, 1) >= 11 && kchg(kc, 2) == 0) {
          goto statement_110;
        }
      }
      //C
      //C...Pack remaining partons/particles. Origin no longer known.
      i1++;
      FEM_DO_SAFE(j, 1, 5) {
        k(i1, j) = k(i, j);
        p(i1, j) = p(i, j);
        v(i1, j) = v(i, j);
      }
      k(i1, 3) = 0;
      statement_110:;
    }
    n = i1;
    //C
    //C...Selective removal of class of entries. New position of retained.
  }
  else if (medit >= 11 && medit <= 15) {
    i1 = 0;
    FEM_DO_SAFE(i, 1, n) {
      k(i, 3) = fem::mod(k(i, 3), mstu(5));
      if (medit == 11 && k(i, 1) < 0) {
        goto statement_120;
      }
      if (medit == 12 && k(i, 1) == 0) {
        goto statement_120;
      }
      if (medit == 13 && (k(i, 1) == 11 || k(i, 1) == 12 || k(i,
          1) == 15) && k(i, 2) != 94) {
        goto statement_120;
      }
      if (medit == 14 && (k(i, 1) == 13 || k(i, 1) == 14 || k(i, 2) == 94)) {
        goto statement_120;
      }
      if (medit == 15 && k(i, 1) >= 21) {
        goto statement_120;
      }
      i1++;
      k(i, 3) += mstu(5) * i1;
      statement_120:;
    }
    //C
    //C...Find new event history information and replace old.
    FEM_DO_SAFE(i, 1, n) {
      if (k(i, 1) <= 0 || k(i, 1) > 20 || k(i, 3) / mstu(5) == 0) {
        goto statement_140;
      }
      id = i;
      statement_130:
      im = fem::mod(k(id, 3), mstu(5));
      if (medit == 13 && im > 0 && im <= n) {
        if ((k(im, 1) == 11 || k(im, 1) == 12 || k(im, 1) == 15) && k(im,
            2) != 94) {
          id = im;
          goto statement_130;
        }
      }
      else if (medit == 14 && im > 0 && im <= n) {
        if (k(im, 1) == 13 || k(im, 1) == 14 || k(im, 2) == 94) {
          id = im;
          goto statement_130;
        }
      }
      k(i, 3) = mstu(5) * (k(i, 3) / mstu(5));
      if (im != 0) {
        k(i, 3) += k(im, 3) / mstu(5);
      }
      if (k(i, 1) != 3 && k(i, 1) != 13 && k(i, 1) != 14) {
        if (k(i, 4) > 0 && k(i, 4) <= mstu(4)) {
          k(i, 4) = k(k(i, 4), 3) / mstu(5);
        }
        if (k(i, 5) > 0 && k(i, 5) <= mstu(4)) {
          k(i, 5) = k(k(i, 5), 3) / mstu(5);
        }
      }
      else {
        kcm = fem::mod(k(i, 4) / mstu(5), mstu(5));
        if (kcm > 0 && kcm <= mstu(4)) {
          kcm = k(kcm, 3) / mstu(5);
        }
        kcd = fem::mod(k(i, 4), mstu(5));
        if (kcd > 0 && kcd <= mstu(4)) {
          kcd = k(kcd, 3) / mstu(5);
        }
        k(i, 4) = fem::pow2(mstu(5)) * (k(i, 4) / fem::pow2(mstu(
          5))) + mstu(5) * kcm + kcd;
        kcm = fem::mod(k(i, 5) / mstu(5), mstu(5));
        if (kcm > 0 && kcm <= mstu(4)) {
          kcm = k(kcm, 3) / mstu(5);
        }
        kcd = fem::mod(k(i, 5), mstu(5));
        if (kcd > 0 && kcd <= mstu(4)) {
          kcd = k(kcd, 3) / mstu(5);
        }
        k(i, 5) = fem::pow2(mstu(5)) * (k(i, 5) / fem::pow2(mstu(
          5))) + mstu(5) * kcm + kcd;
      }
      statement_140:;
    }
    //C
    //C...Pack remaining entries.
    i1 = 0;
    FEM_DO_SAFE(i, 1, n) {
      if (k(i, 3) / mstu(5) == 0) {
        goto statement_160;
      }
      i1++;
      FEM_DO_SAFE(j, 1, 5) {
        k(i1, j) = k(i, j);
        p(i1, j) = p(i, j);
        v(i1, j) = v(i, j);
      }
      k(i1, 3) = fem::mod(k(i1, 3), mstu(5));
      statement_160:;
    }
    n = i1;
    //C
    //C...Save top entries at bottom of LUJETS commonblock.
  }
  else if (medit == 21) {
    if (2 * n >= mstu(4)) {
      luerrm(cmn, 11, "(LUEDIT:) no more memory left in LUJETS");
      return;
    }
    FEM_DO_SAFE(i, 1, n) {
      FEM_DO_SAFE(j, 1, 5) {
        k(mstu(4) - i, j) = k(i, j);
        p(mstu(4) - i, j) = p(i, j);
        v(mstu(4) - i, j) = v(i, j);
      }
    }
    mstu(32) = n;
    //C
    //C...Restore bottom entries of commonblock LUJETS to top.
  }
  else if (medit == 22) {
    FEM_DO_SAFE(i, 1, mstu(32)) {
      FEM_DO_SAFE(j, 1, 5) {
        k(i, j) = k(mstu(4) - i, j);
        p(i, j) = p(mstu(4) - i, j);
        v(i, j) = v(mstu(4) - i, j);
      }
    }
    n = mstu(32);
    //C
    //C...Mark primary entries at top of commonblock LUJETS as untreated.
  }
  else if (medit == 23) {
    i1 = 0;
    FEM_DO_SAFE(i, 1, n) {
      kh = k(i, 3);
      if (kh >= 1) {
        if (k(kh, 1) > 20) {
          kh = 0;
        }
      }
      if (kh != 0) {
        goto statement_200;
      }
      i1++;
      if (k(i, 1) > 10 && k(i, 1) <= 20) {
        k(i, 1) = k(i, 1) - 10;
      }
    }
    statement_200:
    n = i1;
    //C
    //C...Place largest axis along z axis and second largest in xy plane.
  }
  else if (medit == 31 || medit == 32) {
    ludbrb(1, n + mstu(3), 0.f, -ulangl(cmn, p(mstu(61), 1), p(mstu(61),
      2)), 0e0, 0e0, 0e0);
    ludbrb(1, n + mstu(3), -ulangl(cmn, p(mstu(61), 3), p(mstu(61),
      1)), 0.f, 0e0, 0e0, 0e0);
    ludbrb(1, n + mstu(3), 0.f, -ulangl(cmn, p(mstu(61) + 1, 1), p(
      mstu(61) + 1, 2)), 0e0, 0e0, 0e0);
    if (medit == 31) {
      return;
    }
    //C
    //C...Rotate to put slim jet along +z axis.
    FEM_DO_SAFE(is, 1, 2) {
      ns(is) = 0;
      pts(is) = 0.f;
      pls(is) = 0.f;
    }
    FEM_DO_SAFE(i, 1, n) {
      if (k(i, 1) <= 0 || k(i, 1) > 10) {
        goto statement_220;
      }
      if (mstu(41) >= 2) {
        kc = lucomp(cmn, k(i, 2));
        if (kc == 0 || kc == 12 || kc == 14 || kc == 16 || kc == 18) {
          goto statement_220;
        }
        if (mstu(41) >= 3 && kchg(kc, 2) == 0 && luchge(cmn, k(i, 2)) == 0) {
          goto statement_220;
        }
      }
      is = fem::fint(2.f - fem::sign(0.5f, p(i, 3)));
      ns(is)++;
      pts(is) += fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)));
      statement_220:;
    }
    if (ns(1) * fem::pow2(pts(2)) < ns(2) * fem::pow2(pts(1))) {
      ludbrb(1, n + mstu(3), paru(1), 0.f, 0e0, 0e0, 0e0);
    }
    //C
    //C...Rotate to put second largest jet into -z,+x quadrant.
    FEM_DO_SAFE(i, 1, n) {
      if (p(i, 3) >= 0.f) {
        goto statement_230;
      }
      if (k(i, 1) <= 0 || k(i, 1) > 10) {
        goto statement_230;
      }
      if (mstu(41) >= 2) {
        kc = lucomp(cmn, k(i, 2));
        if (kc == 0 || kc == 12 || kc == 14 || kc == 16 || kc == 18) {
          goto statement_230;
        }
        if (mstu(41) >= 3 && kchg(kc, 2) == 0 && luchge(cmn, k(i, 2)) == 0) {
          goto statement_230;
        }
      }
      is = fem::fint(2.f - fem::sign(0.5f, p(i, 1)));
      pls(is) = pls(is) - p(i, 3);
      statement_230:;
    }
    if (pls(2) > pls(1)) {
      ludbrb(1, n + mstu(3), 0.f, paru(1), 0e0, 0e0, 0e0);
    }
  }
  //C
}

struct blockdata_ludata_save
{
  fem::variant_bindings ludatr_bindings;
};

//C
//C*********************************************************************
//C
void
blockdata_ludata(
  common& cmn)
{
  FEM_CMN_SVE(blockdata_ludata);
  // COMMON ludat1
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<float> paru(cmn.paru, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_ref<float> parj(cmn.parj, dimension(200));
  // COMMON ludat2
  arr_ref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_ref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_ref<float> parf(cmn.parf, dimension(2000));
  arr_ref<float, 2> vckm(cmn.vckm, dimension(4, 4));
  // COMMON ludat3
  arr_ref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_ref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_ref<float> brat(cmn.brat, dimension(2000));
  arr_ref<int, 2> kfdp(cmn.kfdp, dimension(2000, 5));
  // COMMON ludat4
  str_arr_ref<1> chaf(cmn.chaf, dimension(500));
  //
  common_variant ludatr(cmn.common_ludatr, sve.ludatr_bindings);
  int i = fem::int0;
  int j = fem::int0;
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<int> mrlu(dimension(6));
      mbr<float> rrlu(dimension(100));
      ludatr.allocate(), mrlu, rrlu;
    }
  }
  arr_ref<int> mrlu(ludatr.bind<int>(), dimension(6));
  /* arr_cref<float> rrlu( */ ludatr.bind<float>() /* , dimension(100)) */ ;
  if (is_called_first_time) {
    {
      fem::data_values data;
      data.values, 0, 0, 0, 9000, 10000, 500, 2000, 0;
      data.values, 0, 2, 6, 1, 1, 0, 1, 1;
      data.values, 0, 0, 0, 0, 2, 10, 0, 0;
      data.values, 1, 10, 0, 0, 0, 0, 0, 0;
      data.values, 0, 0, 0, 0, 0, 0, 0, 0;
      data.values, 2, 2, 1, 4, 2, 1, 1, 0;
      data.values, 0, 0, 25, 24, 0, 1, 0, 0;
      data.values, 0, 0, 0, 0, 0, 0, 0, 0;
      data.values, 0, 0, 0, 0, 0, 0, 40*datum(0), 1;
      data.values, 5, 3, 5, 0, 0, 0, 0, 0;
      data.values, 0, 60*datum(0), 7, 2, 1989, 11, 25, 0;
      data.values, 0, 0, 0, 0, 0, 0, 0, 0;
      data.values, 0, 0, 0, 0, 0, 0;
      data, mstu;
    }
    {
      fem::data_values data;
      data.values, 3.1415927f, 6.2831854f, 0.1973f, 5.068f, 0.3894f,
        2.568f, 4*datum(0.f), 0.001f;
      data.values, 0.09f, 0.01f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 2.0f, 1.0f, 0.25f;
      data.values, 2.5f, 0.05f, 0.f, 0.f, 0.0001f, 0.f, 0.f, 2.5f;
      data.values, 1.5f, 7.0f, 1.0f, 0.5f, 2.0f, 3.2f, 0.f, 0.f;
      data.values, 0.f, 40*datum(0.f), 0.0072974f, 0.230f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.20f, 0.25f, 1.0f, 4.0f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.0f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 70*datum(0.f);
      data, paru;
    }
    {
      fem::data_values data;
      data.values, 1, 3, 0, 0, 0, 0, 0, 0;
      data.values, 0, 0, 1, 2, 0, 1, 0, 0;
      data.values, 0, 0, 0, 0, 2, 1, 1, 2;
      data.values, 1, 0, 0, 0, 0, 0, 0, 0;
      data.values, 0, 0, 0, 0, 0, 0, 0, 0;
      data.values, 1, 2, 4, 2, 5, 0, 1, 0;
      data.values, 0, 0, 0, 3, 0, 0, 0, 0;
      data.values, 0, 0, 0, 0, 40*datum(0), 5, 2, 7;
      data.values, 5, 1, 1, 0, 2, 0, 1, 0;
      data.values, 0, 0, 0, 1, 1, 0, 0, 0;
      data.values, 0, 80*datum(0);
      data, mstj;
    }
    {
      fem::data_values data;
      data.values, 0.10f, 0.30f, 0.40f, 0.05f, 0.50f, 0.50f, 0.50f, 0.f;
      data.values, 0.f, 0.f, 0.50f, 0.60f, 0.75f, 0.f, 0.f, 0.f;
      data.values, 0.f, 1.0f, 1.0f, 0.f, 0.35f, 1.0f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.10f, 1.0f;
      data.values, 0.8f, 1.5f, 0.8f, 2.0f, 0.2f, 2.5f, 0.6f, 2.5f;
      data.values, 0.5f, 0.9f, 0.5f, 0.9f, 0.5f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.77f, 0.77f, 0.77f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 1.0f, 0.f, 4.5f, 0.7f, 0.f, 0.003f;
      data.values, 0.5f, 0.5f, 0.f, 0.f, 0.f, 0.f, 10.f, 1000.f;
      data.values, 100.f, 1000.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.4f, 1.0f, 1.0f, 0.f, 10.f, 10.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.02f, 1.0f, 0.2f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 1.5f, 0.5f, 91.2f, 2.40f, 0.02f, 2.0f, 1.0f, 0.25f;
      data.values, 0.002f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.01f, 0.99f;
      data.values, 0.f, 0.f, 0.2f, 0.f, 60*datum(0.f);
      data, parj;
    }
    {
      fem::data_values data;
      data.values, -1, 2, -1, 2, -1, 2, -1, 2;
      data.values, 2*datum(0), -3, 0, -3, 0, -3, 0, -3;
      data.values, 6*datum(0), 3, 9*datum(0), 3, 2*datum(0), 3, 46*datum(0), 2;
      data.values, -1, 2, -1, 2, 3, 11*datum(0), 3, 0;
      data.values, 2*datum(3), 0, 3, 0, 3, 12*datum(0), 3, 0;
      data.values, 2*datum(3), 0, 3, 0, 3, 12*datum(0), 3, 0;
      data.values, 2*datum(3), 0, 3, 0, 3, 12*datum(0), 3, 0;
      data.values, 2*datum(3), 0, 3, 0, 3, 12*datum(0), 3, 0;
      data.values, 2*datum(3), 0, 3, 0, 3, 12*datum(0), 3, 0;
      data.values, 2*datum(3), 0, 3, 0, 3, 72*datum(0), 3, 0;
      data.values, 3, 28*datum(0), 3, 2*datum(0), 3, 8*datum(0), -3, 8*datum(0);
      data.values, 3, 0, -3, 0, 3, -3, 3*datum(0), 3;
      data.values, 6, 0, 3, 5*datum(0), -3, 0, 3, -3;
      data.values, 0, -3, 4*datum(0), -3, 0, 3, 6, -3;
      data.values, 0, 3, -3, 0, -3, 0, 3, 6;
      data.values, 0, 3, 5*datum(0), -3, 0, 3, -3, 0;
      data.values, -3, 114*datum(0);
      FEM_DO_SAFE(i, 1, 500) {
        data, kchg(i, 1);
      }
    }
    {
      fem::data_values data((values, 8*datum(1), 12*datum(0), 2,
        68*datum(0), -1, 410*datum(0)));
      FEM_DO_SAFE(i, 1, 500) {
        data, kchg(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 8*datum(1), 2*datum(0), 8*datum(1), 5*datum(0), 1,
        9*datum(0), 1, 2*datum(0);
      data.values, 1, 2*datum(0), 1, 41*datum(0), 1, 0, 7*datum(1), 10*datum(0);
      data.values, 9*datum(1), 11*datum(0), 9*datum(1), 11*datum(0),
        9*datum(1), 11*datum(0), 9*datum(1), 11*datum(0);
      data.values, 9*datum(1), 11*datum(0), 9*datum(1), 71*datum(0),
        3*datum(1), 22*datum(0), 1, 5*datum(0);
      data.values, 1, 0, 2*datum(1), 6*datum(0), 1, 0, 2*datum(1), 6*datum(0);
      data.values, 2*datum(1), 0, 5*datum(1), 0, 6*datum(1), 4*datum(0),
        6*datum(1), 4*datum(0);
      data.values, 16*datum(1), 4*datum(0), 6*datum(1), 114*datum(0);
      FEM_DO_SAFE(i, 1, 500) {
        data, kchg(i, 3);
      }
    }
    {
      fem::data_values data;
      data.values, .0099f, .0056f, .199f, 1.35f, 5.f, 90.f, 120.f, 200.f;
      data.values, 2*datum(0.f), .00051f, 0.f, .1057f, 0.f, 1.7841f, 0.f, 60.f;
      data.values, 5*datum(0.f), 91.2f, 80.f, 15.f, 6*datum(0.f),
        300.f, 900.f, 600.f;
      data.values, 300.f, 900.f, 300.f, 2*datum(0.f), 5000.f,
        60*datum(0.f), .1396f, .4977f;
      data.values, .4936f, 1.8693f, 1.8645f, 1.9693f, 5.2794f,
        5.2776f, 5.47972f, 0.f;
      data.values, .135f, .5488f, .9575f, 2.9796f, 9.4f, 117.99f, 238.f, 397.f;
      data.values, 2*datum(0.f), .7669f, .8962f, .8921f, 2.0101f,
        2.0071f, 2.1127f, 2*datum(5.3354f);
      data.values, 5.5068f, 0.f, .77f, .782f, 1.0194f, 3.0969f, 9.4603f, 118.f;
      data.values, 238.f, 397.f, 2*datum(0.f), 1.233f, 2*datum(1.3f),
        2*datum(2.322f), 2.51f, 2*datum(5.73f);
      data.values, 5.97f, 0.f, 1.233f, 1.17f, 1.41f, 3.46f, 9.875f, 118.42f;
      data.values, 238.42f, 397.42f, 2*datum(0.f), .983f, 2*datum(1.429f),
        2*datum(2.272f), 2.46f, 2*datum(5.68f);
      data.values, 5.92f, 0.f, .983f, 1.f, 1.4f, 3.4151f, 9.8598f, 118.4f;
      data.values, 238.4f, 397.4f, 2*datum(0.f), 1.26f, 2*datum(1.401f),
        2*datum(2.372f), 2.56f, 2*datum(5.78f);
      data.values, 6.02f, 0.f, 1.26f, 1.283f, 1.422f, 3.5106f, 9.8919f, 118.5f;
      data.values, 238.5f, 397.5f, 2*datum(0.f), 1.318f, 2*datum(1.426f),
        2*datum(2.422f), 2.61f, 2*datum(5.83f);
      data.values, 6.07f, 0.f, 1.318f, 1.274f, 1.525f, 3.5563f,
        9.9132f, 118.45f;
      data.values, 238.45f, 397.45f, 2*datum(0.f), 2*datum(.4977f),
        83*datum(0.f), 1.1156f, 5*datum(0.f), 2.2849f;
      data.values, 0.f, 2*datum(2.46f), 6*datum(0.f), 5.62f, 0.f,
        2*datum(5.84f), 6*datum(0.f), .9396f;
      data.values, .9383f, 0.f, 1.1974f, 1.1926f, 1.1894f, 1.3213f,
        1.3149f, 0.f;
      data.values, 2.454f, 2.4529f, 2.4522f, 2*datum(2.55f), 2.73f,
        4*datum(0.f), 3*datum(5.8f), 2*datum(5.96f);
      data.values, 6.12f, 4*datum(0.f), 1.234f, 1.233f, 1.232f,
        1.231f, 1.3872f, 1.3837f;
      data.values, 1.3828f, 1.535f, 1.5318f, 1.6724f, 3*datum(2.5f),
        2*datum(2.63f), 2.8f, 4*datum(0.f);
      data.values, 3*datum(5.81f), 2*datum(5.97f), 6.13f, 114*datum(0.f);
      FEM_DO_SAFE(i, 1, 500) {
        data, pmas(i, 1);
      }
    }
    {
      fem::data_values data;
      data.values, 22*datum(0.f), 2.4f, 2.3f, 88*datum(0.f), .0002f,
        .001f, 6*datum(0.f), .149f;
      data.values, .0505f, .0513f, 7*datum(0.f), .153f, .0085f,
        .0044f, 7*datum(0.f), .15f;
      data.values, 2*datum(.09f), 2*datum(.06f), .04f, 3*datum(.1f),
        0.f, .15f, .335f, .08f;
      data.values, 2*datum(.01f), 5*datum(0.f), .057f, 2*datum(.287f),
        2*datum(.06f), .04f, 3*datum(.1f), 0.f;
      data.values, .057f, 0.f, .25f, .0135f, 6*datum(0.f), .4f,
        2*datum(.184f), 2*datum(.06f);
      data.values, .04f, 3*datum(.1f), 0.f, .4f, .025f, .055f,
        .0135f, 6*datum(0.f);
      data.values, .11f, .115f, .099f, 2*datum(.06f), 4*datum(.1f),
        0.f, .11f, .185f;
      data.values, .076f, .0026f, 146*datum(0.f), 4*datum(.115f),
        .039f, 2*datum(.036f), .0099f, .0091f;
      data.values, 131*datum(0.f);
      FEM_DO_SAFE(i, 1, 500) {
        data, pmas(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 22*datum(0.f), 2*datum(20.f), 88*datum(0.f),
        .002f, .005f, 6*datum(0.f), .4f, 2*datum(.2f);
      data.values, 7*datum(0.f), .4f, .1f, .015f, 7*datum(0.f), .25f,
        2*datum(.01f), 3*datum(.08f);
      data.values, 2*datum(.2f), .12f, 0.f, .25f, .2f, .001f, 2*datum(.02f),
        5*datum(0.f);
      data.values, .05f, 2*datum(.4f), 3*datum(.08f), 2*datum(.2f),
        .12f, 0.f, .05f, 0.f;
      data.values, .35f, .05f, 6*datum(0.f), 3*datum(.3f), 2*datum(.08f),
        .06f, 2*datum(.2f), .12f;
      data.values, 0.f, .3f, .05f, .025f, .001f, 6*datum(0.f), .25f,
        4*datum(.12f);
      data.values, 4*datum(.2f), 0.f, .25f, .17f, .2f, .01f,
        146*datum(0.f), 4*datum(.14f);
      data.values, .04f, 2*datum(.035f), 2*datum(.05f), 131*datum(0.f);
      FEM_DO_SAFE(i, 1, 500) {
        data, pmas(i, 3);
      }
    }
    {
      fem::data_values data;
      data.values, 12*datum(0.f), 658650.f, 0.f, .091f, 68*datum(0.f),
        .1f, .43f, 15*datum(0.f);
      data.values, 7803.f, 0.f, 3709.f, .32f, .128f, .131f, 3*datum(.393f),
        84*datum(0.f);
      data.values, .004f, 26*datum(0.f), 15540.f, 26.75f, 83*datum(0.f),
        78.88f, 5*datum(0.f), .054f;
      data.values, 0.f, 2*datum(.13f), 6*datum(0.f), .393f, 0.f,
        2*datum(.393f), 9*datum(0.f), 44.3f;
      data.values, 0.f, 24.f, 49.1f, 86.9f, 6*datum(0.f), .13f,
        9*datum(0.f), .393f;
      data.values, 13*datum(0.f), 24.6f, 130*datum(0.f);
      FEM_DO_SAFE(i, 1, 500) {
        data, pmas(i, 4);
      }
    }
    {
      fem::data_values data;
      data.values, 0.5f, 0.25f, 0.5f, 0.25f, 1.f, 0.5f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.5f, 0.f, 0.5f, 0.f, 1.f, 1.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.5f, 0.f, 0.5f, 0.f;
      data.values, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.5f, 0.f;
      data.values, 0.5f, 0.f, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.5f, 0.f, 0.5f, 0.f, 1.f, 1.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.5f, 0.f, 0.5f, 0.f, 1.f, 1.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.75f, 0.5f, 0.f, 0.1667f;
      data.values, 0.0833f, 0.1667f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 1.f, 0.3333f, 0.6667f, 0.3333f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.325f, 0.325f, 0.5f, 1.6f;
      data.values, 5.0f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.11f;
      data.values, 0.16f, 0.048f, 0.50f, 0.45f, 0.55f, 0.60f, 0.f, 0.f;
      data.values, 0.2f, 0.1f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 0.f, 0.f, 1870*datum(0.f);
      data, parf;
    }
    {
      static const float values[] = {
        0.95150f, 0.04847f, 0.00003f, 0.00000f, 0.04847f, 0.94936f,
          0.00217f, 0.00000f, 0.00003f, 0.00217f, 0.99780f, 0.00000f,
          0.00000f, 0.00000f, 0.00000f, 1.00000f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 4) {
        FEM_DO_SAFE(j, 1, 4) {
          data, vckm(i, j);
        }
      }
    }
    {
      fem::data_values data;
      data.values, 14*datum(0), 1, 0, 1, 5*datum(0), 3*datum(1), 6*datum(0), 1;
      data.values, 4*datum(0), 1, 2*datum(0), 1, 42*datum(0), 7*datum(1),
        12*datum(0), 1;
      data.values, 0, 6*datum(1), 0, 8*datum(1), 2*datum(0), 9*datum(1),
        0, 8*datum(1);
      data.values, 2*datum(0), 9*datum(1), 0, 8*datum(1), 2*datum(0),
        9*datum(1), 0, 8*datum(1);
      data.values, 2*datum(0), 9*datum(1), 0, 8*datum(1), 2*datum(0),
        9*datum(1), 0, 8*datum(1);
      data.values, 3*datum(0), 1, 83*datum(0), 1, 5*datum(0), 1, 0, 2*datum(1);
      data.values, 6*datum(0), 1, 0, 2*datum(1), 9*datum(0), 5*datum(1),
        0, 6*datum(1);
      data.values, 4*datum(0), 6*datum(1), 4*datum(0), 16*datum(1),
        4*datum(0), 6*datum(1), 114*datum(0);
      FEM_DO_SAFE(i, 1, 500) {
        data, mdcy(i, 1);
      }
    }
    {
      fem::data_values data;
      data.values, 1, 9, 17, 25, 33, 41, 49, 57;
      data.values, 2*datum(0), 65, 69, 71, 76, 78, 118, 120;
      data.values, 125, 2*datum(0), 127, 136, 149, 166, 186, 6*datum(0);
      data.values, 203, 4*datum(0), 219, 2*datum(0), 227, 42*datum(0), 236, 237;
      data.values, 241, 250, 252, 254, 256, 11*datum(0), 276, 277;
      data.values, 279, 285, 406, 574, 606, 607, 608, 0;
      data.values, 609, 611, 617, 623, 624, 625, 626, 627;
      data.values, 2*datum(0), 628, 629, 632, 635, 638, 640, 641;
      data.values, 642, 643, 0, 644, 645, 650, 658, 661;
      data.values, 670, 685, 686, 2*datum(0), 687, 688, 693, 698;
      data.values, 700, 702, 703, 705, 707, 0, 709, 710;
      data.values, 713, 717, 718, 719, 721, 722, 2*datum(0), 723;
      data.values, 726, 728, 730, 734, 738, 740, 744, 748;
      data.values, 0, 752, 755, 759, 763, 765, 767, 769;
      data.values, 770, 2*datum(0), 771, 773, 775, 777, 779, 781;
      data.values, 784, 786, 788, 0, 791, 793, 806, 810;
      data.values, 812, 814, 816, 817, 2*datum(0), 818, 824, 835;
      data.values, 846, 854, 862, 867, 875, 883, 0, 888;
      data.values, 895, 903, 905, 907, 909, 911, 912, 2*datum(0);
      data.values, 913, 921, 83*datum(0), 923, 5*datum(0), 927, 0, 1001;
      data.values, 1002, 6*datum(0), 1003, 0, 1004, 1005, 9*datum(0), 1006;
      data.values, 1008, 1009, 1012, 1013, 0, 1015, 1016, 1017;
      data.values, 1018, 1019, 1020, 4*datum(0), 1021, 1022, 1023, 1024;
      data.values, 1025, 1026, 4*datum(0), 1027, 1028, 1031, 1034, 1035;
      data.values, 1038, 1041, 1044, 1046, 1048, 1052, 1053, 1054;
      data.values, 1055, 1057, 1059, 4*datum(0), 1060, 1061, 1062, 1063;
      data.values, 1064, 1065, 114*datum(0);
      FEM_DO_SAFE(i, 1, 500) {
        data, mdcy(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 8*datum(8), 2*datum(0), 4, 2, 5, 2, 40, 2;
      data.values, 5, 2, 2*datum(0), 9, 13, 17, 20, 17;
      data.values, 6*datum(0), 16, 4*datum(0), 8, 2*datum(0), 9, 42*datum(0), 1;
      data.values, 4, 9, 3*datum(2), 20, 11*datum(0), 1, 2, 6;
      data.values, 121, 168, 32, 3*datum(1), 0, 2, 2*datum(6), 5*datum(1);
      data.values, 2*datum(0), 1, 3*datum(3), 2, 4*datum(1), 0, 1, 5;
      data.values, 8, 3, 9, 15, 2*datum(1), 2*datum(0), 1, 2*datum(5);
      data.values, 2*datum(2), 1, 3*datum(2), 0, 1, 3, 4, 2*datum(1);
      data.values, 2, 2*datum(1), 2*datum(0), 3, 2*datum(2), 2*datum(4),
        2, 3*datum(4);
      data.values, 0, 3, 2*datum(4), 3*datum(2), 2*datum(1), 2*datum(0),
        5*datum(2), 3;
      data.values, 2*datum(2), 3, 0, 2, 13, 4, 3*datum(2), 2*datum(1);
      data.values, 2*datum(0), 6, 2*datum(11), 2*datum(8), 5, 2*datum(8), 5, 0;
      data.values, 7, 8, 4*datum(2), 2*datum(1), 2*datum(0), 8, 2, 83*datum(0);
      data.values, 4, 5*datum(0), 74, 0, 2*datum(1), 6*datum(0), 1, 0;
      data.values, 2*datum(1), 9*datum(0), 2, 1, 3, 1, 2, 0;
      data.values, 6*datum(1), 4*datum(0), 6*datum(1), 4*datum(0), 1,
        2*datum(3), 1, 3*datum(3);
      data.values, 2*datum(2), 4, 3*datum(1), 2*datum(2), 1, 4*datum(0),
        6*datum(1), 114*datum(0);
      FEM_DO_SAFE(i, 1, 500) {
        data, mdcy(i, 3);
      }
    }
    {
      fem::data_values data;
      data.values, 6*datum(1), -1, 7*datum(1), -1, 7*datum(1), -1,
        7*datum(1), -1;
      data.values, 7*datum(1), -1, 7*datum(1), -1, 85*datum(1),
        2*datum(-1), 7*datum(1), 2*datum(-1);
      data.values, 3*datum(1), 2*datum(-1), 6*datum(1), 2*datum(-1),
        6*datum(1), 3*datum(-1), 3*datum(1), -1;
      data.values, 3*datum(1), -1, 3*datum(1), 5*datum(-1), 3*datum(1),
        -1, 6*datum(1), 2*datum(-1);
      data.values, 3*datum(1), -1, 11*datum(1), 2*datum(-1), 6*datum(1),
        2*datum(-1), 3*datum(1), -1;
      data.values, 3*datum(1), -1, 4*datum(1), 2*datum(-1), 2*datum(1),
        -1, 488*datum(1), 2*datum(0);
      data.values, 1275*datum(1);
      FEM_DO_SAFE(i, 1, 2000) {
        data, mdme(i, 1);
      }
    }
    {
      fem::data_values data;
      data.values, 70*datum(102), 42, 6*datum(102), 2*datum(42),
        2*datum(0), 7*datum(41), 2*datum(0), 23*datum(41);
      data.values, 6*datum(102), 45, 28*datum(102), 8*datum(32),
        9*datum(0), 16*datum(32), 4*datum(0), 8*datum(32);
      data.values, 4*datum(0), 32, 4*datum(0), 8*datum(32), 8*datum(0),
        4*datum(32), 4*datum(0), 6*datum(32);
      data.values, 3*datum(0), 12, 2*datum(42), 2*datum(11), 9*datum(42),
        6*datum(45), 20*datum(46), 7*datum(0);
      data.values, 34*datum(42), 86*datum(0), 2*datum(25), 26,
        24*datum(42), 142*datum(0), 25, 26;
      data.values, 0, 10*datum(42), 19*datum(0), 2*datum(13), 3*datum(85),
        0, 2, 4*datum(0);
      data.values, 2, 8*datum(0), 2*datum(32), 87, 88, 3*datum(3), 0,
        2*datum(3);
      data.values, 0, 2*datum(3), 0, 3, 5*datum(0), 3, 1, 0;
      data.values, 3, 2*datum(0), 2*datum(3), 3*datum(0), 1, 4*datum(0),
        12, 3*datum(0);
      data.values, 4*datum(32), 2*datum(4), 6*datum(0), 5*datum(32),
        2*datum(4), 2*datum(45), 87, 88;
      data.values, 30*datum(0), 12, 32, 0, 32, 87, 88, 41*datum(0);
      data.values, 12, 0, 32, 0, 32, 87, 88, 40*datum(0);
      data.values, 12, 0, 32, 0, 32, 87, 88, 88*datum(0);
      data.values, 12, 0, 32, 0, 32, 87, 88, 2*datum(0);
      data.values, 4*datum(42), 8*datum(0), 14*datum(42), 50*datum(0),
        10*datum(13), 2*datum(84), 3*datum(85), 14*datum(0);
      data.values, 84, 5*datum(0), 85, 974*datum(0);
      FEM_DO_SAFE(i, 1, 2000) {
        data, mdme(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 70*datum(0.f), 1.f, 6*datum(0.f), 2*datum(.177f),
        .108f, .225f, .003f, .06f;
      data.values, .02f, .025f, .013f, 2*datum(.004f), .007f, .014f,
        2*datum(.002f), 2*datum(.001f);
      data.values, .054f, .014f, .016f, .005f, 2*datum(.012f),
        5*datum(.006f), .002f, 2*datum(.001f);
      data.values, 5*datum(.002f), 6*datum(0.f), 1.f, 28*datum(0.f),
        .143f, .111f, .143f, .111f;
      data.values, .143f, .085f, 2*datum(0.f), .03f, .058f, .03f, .058f, .03f;
      data.values, .058f, 3*datum(0.f), .25f, .01f, 2*datum(0.f),
        .01f, .25f, 4*datum(0.f);
      data.values, .24f, 5*datum(0.f), 3*datum(.08f), 3*datum(0.f),
        .01f, .08f, .82f, 5*datum(0.f);
      data.values, .09f, 6*datum(0.f), .143f, .111f, .143f, .111f, .143f, .085f;
      data.values, 2*datum(0.f), .03f, .058f, .03f, .058f, .03f,
        .058f, 4*datum(0.f);
      data.values, 1.f, 5*datum(0.f), 4*datum(.215f), 2*datum(0.f),
        2*datum(.07f), 0.f, 1.f, 2*datum(.08f);
      data.values, .76f, .08f, 2*datum(.112f), .05f, .476f, .08f, .14f, .01f;
      data.values, .015f, .005f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f;
      data.values, .25f, .01f, 2*datum(0.f), .01f, .25f, 4*datum(0.f),
        .24f, 5*datum(0.f);
      data.values, 3*datum(.08f), 0.f, 1.f, 2*datum(.5f), .635f,
        .212f, .056f, .017f;
      data.values, .048f, .032f, .035f, .03f, 2*datum(.015f), .044f,
        2*datum(.022f), 9*datum(.001f);
      data.values, .035f, .03f, 2*datum(.015f), .044f, 2*datum(.022f),
        9*datum(.001f), .028f, .017f;
      data.values, .066f, .02f, .008f, 2*datum(.006f), .003f, .001f,
        2*datum(.002f), .003f;
      data.values, .001f, 2*datum(.002f), .005f, .002f, .005f, .006f,
        .004f, .012f;
      data.values, 2*datum(.005f), .008f, 2*datum(.005f), .037f,
        .004f, .067f, 2*datum(.01f), 2*datum(.001f);
      data.values, 3*datum(.002f), .003f, 8*datum(.002f), .005f,
        4*datum(.004f), .015f, .005f, .027f;
      data.values, 2*datum(.005f), .007f, .014f, .007f, .01f, .008f,
        .012f, .015f;
      data.values, 11*datum(.002f), 3*datum(.004f), .002f, .004f,
        6*datum(.002f), 2*datum(.004f), .005f, .011f;
      data.values, .005f, .015f, .02f, 2*datum(.01f), 3*datum(.004f),
        5*datum(.002f), .015f, .02f;
      data.values, 2*datum(.01f), 3*datum(.004f), 5*datum(.002f),
        .038f, .048f, .082f, .06f, .028f;
      data.values, .021f, 2*datum(.005f), 2*datum(.002f), .005f,
        .018f, .005f, .01f, .008f;
      data.values, .005f, 3*datum(.004f), .001f, 3*datum(.003f),
        .001f, 2*datum(.002f), .003f, 2*datum(.002f);
      data.values, 2*datum(.001f), .002f, .001f, .002f, .001f, .005f,
        4*datum(.003f), .001f;
      data.values, 2*datum(.002f), .003f, 2*datum(.001f), .013f,
        .03f, .058f, .055f, 3*datum(.003f);
      data.values, 2*datum(.01f), .007f, .019f, 4*datum(.005f),
        .015f, 3*datum(.005f), 8*datum(.002f), 3*datum(.001f);
      data.values, .002f, 2*datum(.001f), .003f, 16*datum(.001f);
      FEM_DO_SAFE(i, 1, 525) {
        data, brat(i);
      }
    }
    {
      fem::data_values data;
      data.values, .019f, 2*datum(.003f), .002f, .005f, .004f, .008f,
        .003f, .006f;
      data.values, .003f, .01f, 5*datum(.002f), 2*datum(.001f),
        2*datum(.002f), 11*datum(.001f), .002f, 14*datum(.001f);
      data.values, .018f, .005f, .01f, 2*datum(.015f), .017f, 4*datum(.015f),
        .017f, 3*datum(.015f);
      data.values, .025f, .08f, 2*datum(.025f), .04f, .001f, 2*datum(.005f),
        .02f, .04f;
      data.values, 2*datum(.06f), .04f, .01f, 4*datum(.005f), .25f,
        .115f, 3*datum(1.f), .988f;
      data.values, .012f, .389f, .319f, .237f, .049f, .005f, .001f, .441f;
      data.values, .205f, .301f, .03f, .022f, .001f, 6*datum(1.f), .665f, .333f;
      data.values, .002f, .666f, .333f, .001f, .49f, .34f, .17f, .52f;
      data.values, .48f, 5*datum(1.f), .893f, .08f, .017f, 2*datum(.005f),
        .495f, .343f;
      data.values, 3*datum(.043f), .019f, .013f, .001f, 2*datum(.069f),
        .862f, 3*datum(.027f), .015f;
      data.values, .045f, .015f, .045f, .77f, .029f, 6*datum(.02f),
        5*datum(.05f), .115f;
      data.values, .015f, .5f, 0.f, 3*datum(1.f), .28f, .14f, .313f, .157f;
      data.values, .11f, .28f, .14f, .313f, .157f, .11f, .667f, .333f;
      data.values, .667f, .333f, 1.f, .667f, .333f, .667f, .333f, 2*datum(.5f);
      data.values, 1.f, .333f, .334f, .333f, 4*datum(.25f), 2*datum(1.f),
        .3f, .7f;
      data.values, 2*datum(1.f), .8f, 2*datum(.1f), .667f, .333f,
        .667f, .333f, .6f;
      data.values, .3f, .067f, .033f, .6f, .3f, .067f, .033f, 2*datum(.5f);
      data.values, .6f, .3f, .067f, .033f, .6f, .3f, .067f, .033f;
      data.values, 2*datum(.4f), 2*datum(.1f), .8f, 2*datum(.1f),
        .52f, .26f, 2*datum(.11f), .62f;
      data.values, .31f, 2*datum(.035f), .007f, .993f, .02f, .98f, .3f, .7f;
      data.values, 2*datum(1.f), 2*datum(.5f), .667f, .333f, .667f,
        .333f, .667f, .333f;
      data.values, .667f, .333f, 2*datum(.35f), .3f, .667f, .333f, .667f, .333f;
      data.values, 2*datum(.35f), .3f, 2*datum(.5f), 3*datum(.14f),
        .1f, .05f, 4*datum(.08f), .028f;
      data.values, .027f, .028f, .027f, 4*datum(.25f), .273f, .727f, .35f, .65f;
      data.values, .3f, .7f, 2*datum(1.f), 2*datum(.35f), .144f,
        .105f, .048f, .003f;
      data.values, .332f, .166f, .168f, .084f, .086f, .043f, .059f,
        2*datum(.029f);
      data.values, 2*datum(.002f), .332f, .166f, .168f, .084f, .086f,
        .043f, .059f;
      data.values, 2*datum(.029f), 2*datum(.002f), .3f, .15f, .16f,
        .08f, .13f, .06f;
      data.values, .08f, .04f, .3f, .15f, .16f, .08f, .13f, .06f;
      data.values, .08f, .04f, 2*datum(.4f), .1f, 2*datum(.05f), .3f,
        .15f, .16f;
      data.values, .08f, .13f, .06f, .08f, .04f, .3f, .15f, .16f;
      data.values, .08f, .13f, .06f, .08f, .04f, 2*datum(.4f), .1f,
        2*datum(.05f);
      data.values, 2*datum(.35f), .144f, .105f, 2*datum(.024f);
      FEM_DO_SAFE(i, 526, 893) {
        data, brat(i);
      }
    }
    {
      fem::data_values data;
      data.values, .003f, .573f, .287f, .063f, .028f, 2*datum(.021f),
        .004f, .003f;
      data.values, 2*datum(.5f), .15f, .85f, .22f, .78f, .3f, .7f, 2*datum(1.f);
      data.values, .217f, .124f, 2*datum(.193f), 2*datum(.135f),
        .002f, .001f, .686f, .314f;
      data.values, .641f, .357f, 2*datum(.001f), .018f, 2*datum(.005f),
        .003f, .002f, 2*datum(.006f);
      data.values, .018f, 2*datum(.005f), .003f, .002f, 2*datum(.006f),
        .005f, .025f, .015f;
      data.values, .006f, 2*datum(.005f), .004f, .005f, 5*datum(.004f),
        2*datum(.002f), 2*datum(.004f), .003f;
      data.values, .002f, 2*datum(.003f), 3*datum(.002f), 2*datum(.001f),
        .002f, 2*datum(.001f), 2*datum(.002f), 5*datum(.001f);
      data.values, 4*datum(.003f), 2*datum(.005f), 2*datum(.002f),
        2*datum(.001f), 2*datum(.002f), 2*datum(.001f), .255f, .057f;
      data.values, 2*datum(.035f), .15f, 2*datum(.075f), .03f,
        2*datum(.015f), 5*datum(1.f), .999f, .001f;
      data.values, 1.f, .516f, .483f, .001f, 1.f, .995f, .005f, 13*datum(1.f);
      data.values, .331f, .663f, .006f, .663f, .331f, .006f, 1.f, .88f;
      data.values, 2*datum(.06f), .88f, 2*datum(.06f), .88f, 2*datum(.06f),
        .667f, 2*datum(.333f), .667f;
      data.values, .676f, .234f, .085f, .005f, 3*datum(1.f), 4*datum(.5f),
        7*datum(1.f), 935*datum(0.f);
      FEM_DO_SAFE(i, 894, 2000) {
        data, brat(i);
      }
    }
    {
      fem::data_values data;
      data.values, 21, 22, 23, 4*datum(-24), 25, 21, 22, 23;
      data.values, 4*datum(24), 25, 21, 22, 23, 4*datum(-24), 25, 21;
      data.values, 22, 23, 4*datum(24), 25, 21, 22, 23, 4*datum(-24);
      data.values, 25, 21, 22, 23, 4*datum(24), 25, 21, 22;
      data.values, 23, 4*datum(-24), 25, 21, 22, 23, 4*datum(24), 25;
      data.values, 22, 23, -24, 25, 23, 24, -12, 22;
      data.values, 23, -24, 25, 23, 24, -12, -14, 34*datum(16);
      data.values, 22, 23, -24, 25, 23, 24, -89, 22;
      data.values, 23, -24, 25, 23, 24, 1, 2, 3;
      data.values, 4, 5, 6, 7, 8, 21, 1, 2;
      data.values, 3, 4, 5, 6, 7, 8, 11, 13;
      data.values, 15, 17, 37, 1, 2, 3, 4, 5;
      data.values, 6, 7, 8, 11, 12, 13, 14, 15;
      data.values, 16, 17, 18, 37, 4*datum(-1), 4*datum(-3), 4*datum(-5),
        4*datum(-7);
      data.values, -11, -13, -15, -17, 1, 2, 3, 4;
      data.values, 5, 6, 7, 8, 11, 13, 15, 17;
      data.values, 21, 2*datum(22), 23, 24, 1, 2, 3, 4;
      data.values, 5, 6, 7, 8, 11, 12, 13, 14;
      data.values, 15, 16, 17, 18, -1, -3, -5, -7;
      data.values, -11, -13, -15, -17, 1, 2, 3, 4;
      data.values, 5, 6, 11, 13, 15, 82, -11, -13;
      data.values, 2*datum(2), -12, -14, -16, 2*datum(-2), 2*datum(-4), -2, -4;
      data.values, 2*datum(89), 2*datum(-89), 2*datum(89), 4*datum(-1),
        4*datum(-3), 4*datum(-5), 4*datum(-7), -11;
      data.values, -13, -15, -17, -13, 130, 310, -13, 3*datum(211);
      data.values, 12, 14, 16*datum(-11), 16*datum(-13), -311, -313, -311, -313;
      data.values, -311, -313, -311, -313, 2*datum(111), 2*datum(221),
        2*datum(331), 2*datum(113);
      data.values, 2*datum(223), 2*datum(333), -311, -313, 2*datum(-311),
        -313, 3*datum(-311), -321;
      data.values, -323, -321, 2*datum(211), 2*datum(213), -213, 113,
        3*datum(213), 3*datum(211);
      data.values, 2*datum(213), 2*datum(-311), -313, -321, 2*datum(-311),
        -313, -311, -313;
      data.values, 4*datum(-311), -321, -323, 2*datum(-321), 3*datum(211),
        213, 2*datum(211), 213;
      data.values, 5*datum(211), 213, 4*datum(211), 3*datum(213),
        211, 213, 321, 311;
      data.values, 3, 2*datum(2), 12*datum(-11), 12*datum(-13), -321,
        -323, -321, -323;
      data.values, -311, -313, -311, -313, -311, -313, -311, -313;
      data.values, -311, -313, -311, -321, -323, -321, -323, 211;
      data.values, 213, 211, 213, 111, 221, 331, 113, 223;
      data.values, 333, 221, 331, 113, 223, 113, 223, 113;
      data.values, 223, 333, 223, 333, 321, 323, 321, 323;
      data.values, 311, 313, -321, -323, 3*datum(-321), -323, 2*datum(-321),
        -323;
      data.values, -321, -311, -313, 3*datum(-311), -313, 2*datum(-311),
        -313, -321;
      data.values, -323, 3*datum(-321);
      FEM_DO_SAFE(i, 1, 499) {
        data, kfdp(i, 1);
      }
    }
    {
      fem::data_values data;
      data.values, -323, 2*datum(-321), -311, 2*datum(333), 211, 213,
        2*datum(211), 2*datum(213);
      data.values, 4*datum(211), 10*datum(111), -321, -323, 5*datum(-321),
        -323, 2*datum(-321), -311;
      data.values, -313, 4*datum(-311), -313, 4*datum(-311), -321,
        -323, 2*datum(-321), -323;
      data.values, -321, -313, -311, -313, -311, 211, 213, 2*datum(211);
      data.values, 213, 4*datum(211), 111, 221, 113, 223, 113, 223;
      data.values, 2*datum(3), -15, 5*datum(-11), 5*datum(-13), 221,
        331, 333, 221;
      data.values, 331, 333, 211, 213, 211, 213, 321, 323;
      data.values, 321, 323, 2212, 221, 331, 333, 221, 2*datum(2);
      data.values, 3*datum(0), 3*datum(22), 111, 211, 2*datum(22),
        2*datum(211), 111, 3*datum(22);
      data.values, 111, 3*datum(21), 2*datum(0), 211, 321, 3*datum(311),
        2*datum(321), 421;
      data.values, 2*datum(411), 2*datum(421), 431, 511, 521, 531,
        2*datum(211), 22;
      data.values, 211, 2*datum(111), 321, 130, -213, 113, 213, 211;
      data.values, 22, 111, 11, 13, 82, 11, 13, 15;
      data.values, 1, 2, 3, 4, 21, 22, 11, 12;
      data.values, 13, 14, 15, 16, 1, 2, 3, 4;
      data.values, 5, 21, 22, 2*datum(89), 2*datum(0), 223, 321, 311;
      data.values, 323, 313, 2*datum(311), 321, 313, 323, 321, 421;
      data.values, 2*datum(411), 421, 433, 521, 2*datum(511), 521, 523, 513;
      data.values, 223, 213, 113, -213, 313, -313, 323, -323;
      data.values, 82, 21, 663, 21, 2*datum(0), 221, 213, 113;
      data.values, 321, 2*datum(311), 321, 421, 411, 423, 413, 411;
      data.values, 421, 413, 423, 431, 433, 521, 511, 523;
      data.values, 513, 511, 521, 513, 523, 521, 511, 531;
      data.values, 533, 221, 213, -213, 211, 111, 321, 130;
      data.values, 211, 111, 321, 130, 443, 82, 553, 21;
      data.values, 663, 21, 2*datum(0), 113, 213, 323, 2*datum(313), 323;
      data.values, 423, 2*datum(413), 423, 421, 411, 433, 523, 2*datum(513);
      data.values, 523, 521, 511, 533, 213, -213, 10211, 10111;
      data.values, -10211, 2*datum(221), 213, 2*datum(113), -213,
        2*datum(321), 2*datum(311), 313;
      data.values, -313, 323, -323, 443, 82, 553, 21, 663;
      data.values, 21, 2*datum(0), 213, 113, 221, 223, 321, 211;
      data.values, 321, 311, 323, 313, 323, 313, 321, 5*datum(311);
      data.values, 321, 313, 323, 313, 323, 311, 4*datum(321), 421;
      data.values, 411, 423, 413, 423, 413, 421, 2*datum(411), 421;
      data.values, 413, 423, 413, 423, 411, 2*datum(421), 411, 433;
      data.values, 2*datum(431), 521, 511, 523, 513, 523, 513, 521;
      FEM_DO_SAFE(i, 500, 873) {
        data, kfdp(i, 1);
      }
    }
    {
      fem::data_values data;
      data.values, 2*datum(511), 521, 513, 523, 513, 523, 511, 2*datum(521);
      data.values, 511, 533, 2*datum(531), 213, -213, 221, 223, 321;
      data.values, 130, 111, 211, 111, 2*datum(211), 321, 130, 221;
      data.values, 111, 321, 130, 443, 82, 553, 21, 663;
      data.values, 21, 2*datum(0), 111, 211, -12, 12, -14, 14;
      data.values, 211, 111, 211, 111, 2212, 2*datum(2112), -12, 7*datum(-11);
      data.values, 7*datum(-13), 2*datum(2224), 2*datum(2212),
        2*datum(2214), 2*datum(3122), 2*datum(3212), 2*datum(3214),
        5*datum(3222);
      data.values, 4*datum(3224), 2*datum(3322), 3324, 2*datum(2224),
        5*datum(2212), 5*datum(2214), 2*datum(2112), 2*datum(2114);
      data.values, 2*datum(3122), 2*datum(3212), 2*datum(3214),
        2*datum(3222), 2*datum(3224), 4*datum(2), 3, 2*datum(2);
      data.values, 1, 2*datum(2), 5*datum(0), 2112, -12, 3122, 2212, 2112;
      data.values, 2212, 3*datum(3122), 3*datum(4122), 4132, 4232, 0,
        3*datum(5122), 5132;
      data.values, 5232, 0, 2112, 2212, 2*datum(2112), 2212, 2112,
        2*datum(2212);
      data.values, 3122, 3212, 3112, 3122, 3222, 3112, 3122, 3222;
      data.values, 3212, 3322, 3312, 3322, 3312, 3122, 3322, 3312;
      data.values, -12, 3*datum(4122), 2*datum(4132), 2*datum(4232),
        4332, 3*datum(5122), 5132, 5232;
      data.values, 5332, 935*datum(0);
      FEM_DO_SAFE(i, 874, 2000) {
        data, kfdp(i, 1);
      }
    }
    {
      fem::data_values data;
      data.values, 3*datum(1), 2, 4, 6, 8, 1, 3*datum(2), 1;
      data.values, 3, 5, 7, 2, 3*datum(3), 2, 4, 6;
      data.values, 8, 3, 3*datum(4), 1, 3, 5, 7, 4;
      data.values, 3*datum(5), 2, 4, 6, 8, 5, 3*datum(6), 1;
      data.values, 3, 5, 7, 6, 3*datum(7), 2, 4, 6;
      data.values, 8, 7, 3*datum(8), 1, 3, 5, 7, 8;
      data.values, 2*datum(11), 12, 11, 12, 2*datum(11), 2*datum(13), 14, 13;
      data.values, 14, 13, 11, 13, -211, -213, -211, -213;
      data.values, -211, -213, 3*datum(-211), -321, -323, -321, -323,
        2*datum(-321);
      data.values, 4*datum(-211), -213, -211, -213, -211, -213, -211, -213;
      data.values, -211, -213, 6*datum(-211), 2*datum(15), 16, 15, 16, 15;
      data.values, 18, 2*datum(17), 18, 17, 18, 17, -1, -2;
      data.values, -3, -4, -5, -6, -7, -8, 21, -1;
      data.values, -2, -3, -4, -5, -6, -7, -8, -11;
      data.values, -13, -15, -17, -37, -1, -2, -3, -4;
      data.values, -5, -6, -7, -8, -11, -12, -13, -14;
      data.values, -15, -16, -17, -18, -37, 2, 4, 6;
      data.values, 8, 2, 4, 6, 8, 2, 4, 6;
      data.values, 8, 2, 4, 6, 8, 12, 14, 16;
      data.values, 18, -1, -2, -3, -4, -5, -6, -7;
      data.values, -8, -11, -13, -15, -17, 21, 22, 2*datum(23);
      data.values, -24, -1, -2, -3, -4, -5, -6, -7;
      data.values, -8, -11, -12, -13, -14, -15, -16, -17;
      data.values, -18, 2, 4, 6, 8, 12, 14, 16;
      data.values, 18, -3, -4, -5, -6, -7, -8, -13;
      data.values, -15, -17, -82, 12, 14, -1, -3, 11;
      data.values, 13, 15, 1, 4, 3, 4, 1, 3;
      data.values, 5, 3, 6, 4, 7, 5, 2, 4;
      data.values, 6, 8, 2, 4, 6, 8, 2, 4;
      data.values, 6, 8, 2, 4, 6, 8, 12, 14;
      data.values, 16, 18, 14, 2*datum(0), 14, 111, 211, 111;
      data.values, -11, -13, 16*datum(12), 16*datum(14), 2*datum(211),
        2*datum(213), 2*datum(321), 2*datum(323);
      data.values, 211, 213, 211, 213, 211, 213, 211, 213;
      data.values, 211, 213, 211, 213, 2*datum(211), 213, 7*datum(211), 213;
      data.values, 211, 111, 211, 111, 2*datum(211), -213, 213, 2*datum(113);
      data.values, 223, 2*datum(113), 221, 321, 2*datum(311), 321,
        313, 4*datum(211);
      data.values, 213, 113, 213, -213, 2*datum(211), 213, 113, 111;
      data.values, 221, 331, 111, 113, 223, 4*datum(113), 223, 6*datum(211);
      data.values, 213, 4*datum(211), -321, -311, 3*datum(-1),
        12*datum(12), 12*datum(14), 2*datum(211);
      data.values, 2*datum(213), 2*datum(111), 2*datum(221), 2*datum(331),
        2*datum(113), 2*datum(223), 333, 2*datum(321);
      data.values, 2*datum(323), 2*datum(-211), 2*datum(-213),
        6*datum(111), 4*datum(221), 2*datum(331), 3*datum(113),
        2*datum(223);
      data.values, 2*datum(-211), 2*datum(-213), 113, 111, 2*datum(211),
        213, 6*datum(211), 321;
      data.values, 2*datum(211), 213, 211, 2*datum(111), 113, 2*datum(223),
        2*datum(321);
      FEM_DO_SAFE(i, 1, 496) {
        data, kfdp(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 323, 321, 2*datum(311), 313, 2*datum(311), 111,
        211, 2*datum(-211);
      data.values, -213, -211, -213, -211, -213, 3*datum(-211),
        5*datum(111), 2*datum(113);
      data.values, 223, 113, 223, 2*datum(211), 213, 5*datum(211),
        213, 3*datum(211);
      data.values, 213, 2*datum(211), 2*datum(111), 221, 113, 223,
        3*datum(321), 323;
      data.values, 2*datum(321), 323, 311, 313, 311, 313, 3*datum(211),
        2*datum(-211);
      data.values, -213, 3*datum(-211), 4*datum(111), 2*datum(113),
        2*datum(-1), 16, 5*datum(12), 5*datum(14);
      data.values, 3*datum(211), 3*datum(213), 2*datum(111), 2*datum(113),
        2*datum(-311), 2*datum(-313), -2112, 3*datum(321);
      data.values, 323, 2*datum(-1), 3*datum(0), 22, 11, 22, 111, -211;
      data.values, 211, 11, 2*datum(-211), 111, 113, 223, 22, 111;
      data.values, 3*datum(21), 2*datum(0), 111, -211, 111, 22, 211, 111;
      data.values, 22, 211, 111, 22, 111, 5*datum(22), 2*datum(-211), 111;
      data.values, -211, 2*datum(111), -321, 310, 211, 111, 2*datum(-211), 221;
      data.values, 22, -11, -13, -82, -11, -13, -15, -1;
      data.values, -2, -3, -4, 2*datum(21), -11, -12, -13, -14;
      data.values, -15, -16, -1, -2, -3, -4, -5, 2*datum(21);
      data.values, 5, 3, 2*datum(0), 211, -213, 113, -211, 111;
      data.values, 223, 211, 111, 211, 111, 223, 211, 111;
      data.values, -211, 2*datum(111), -211, 111, 211, 111, -321, -311;
      data.values, 111, -211, 111, 211, -311, 311, -321, 321;
      data.values, -82, 21, 22, 21, 2*datum(0), 211, 111, 211;
      data.values, -211, 111, 211, 111, 211, 111, 211, 111;
      data.values, -211, 111, -211, 3*datum(111), -211, 111, -211, 111;
      data.values, 211, 111, 211, 111, -321, -311, 3*datum(111), -211;
      data.values, 211, -211, 111, -321, 310, -211, 111, -321;
      data.values, 310, 22, -82, 22, 21, 22, 21, 2*datum(0);
      data.values, 211, 111, -211, 111, 211, 111, 211, 111;
      data.values, -211, 111, 321, 311, 111, -211, 111, 211;
      data.values, 111, -321, -311, 111, -211, 211, -211, 111;
      data.values, 2*datum(211), 111, -211, 211, 111, 211, -321, 2*datum(-311);
      data.values, -321, -311, 311, -321, 321, 22, -82, 22;
      data.values, 21, 22, 21, 2*datum(0), 111, 3*datum(211), -311, 22;
      data.values, -211, 111, -211, 111, -211, 211, -213, 113;
      data.values, 223, 221, 22, 211, 111, 211, 111, 2*datum(211);
      data.values, 213, 113, 223, 221, 22, 211, 111, 211;
      data.values, 111, 4*datum(211), -211, 111, -211, 111, -211, 211;
      data.values, -211, 211, 321, 311;
      FEM_DO_SAFE(i, 497, 863) {
        data, kfdp(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 2*datum(111), 211, -211, 111, -211, 111, -211, 211;
      data.values, -211, 2*datum(211), 111, 211, 111, 4*datum(211), -321, -311;
      data.values, 2*datum(111), 211, -211, 211, 111, 211, -321, 310;
      data.values, 22, -211, 111, 2*datum(-211), -321, 310, 221, 111;
      data.values, -321, 310, 22, -82, 22, 21, 22, 21;
      data.values, 2*datum(0), 111, -211, 11, -11, 13, -13, -211;
      data.values, 111, -211, 111, -211, 111, 22, 11, 7*datum(12);
      data.values, 7*datum(14), -321, -323, -311, -313, -311, -313, 211;
      data.values, 213, 211, 213, 211, 213, 111, 221, 331;
      data.values, 113, 223, 111, 221, 113, 223, 321, 323;
      data.values, 321, -211, -213, 111, 221, 331, 113, 223;
      data.values, 111, 221, 331, 113, 223, 211, 213, 211;
      data.values, 213, 321, 323, 321, 323, 321, 323, 311;
      data.values, 313, 311, 313, 2*datum(-1), -3, -1, 2203, 2*datum(3201);
      data.values, 2203, 2101, 2103, 5*datum(0), -211, 11, 22, 111;
      data.values, 211, 22, -211, 111, 22, -211, 111, 211;
      data.values, 2*datum(22), 0, -211, 111, 211, 2*datum(22), 0,
        2*datum(-211);
      data.values, 111, 22, 111, 211, 22, 211, 2*datum(-211), 2*datum(111);
      data.values, -211, 2*datum(211), 111, 211, -211, 2*datum(111), 211, -321;
      data.values, -211, 111, 11, -211, 111, 211, 111, 22;
      data.values, 111, 2*datum(22), -211, 111, 211, 3*datum(22), 935*datum(0);
      FEM_DO_SAFE(i, 864, 2000) {
        data, kfdp(i, 2);
      }
    }
    {
      fem::data_values data;
      data.values, 70*datum(0), 14, 6*datum(0), 2*datum(16), 2*datum(0),
        5*datum(111), 310, 130;
      data.values, 2*datum(0), 2*datum(111), 310, 130, 113, 211, 223, 221;
      data.values, 2*datum(113), 2*datum(211), 2*datum(223), 2*datum(221),
        2*datum(113), 221, 113, 2*datum(213);
      data.values, -213, 123*datum(0), 4*datum(3), 4*datum(4), 1, 4,
        3, 2*datum(2);
      data.values, 6*datum(81), 25*datum(0), -211, 3*datum(111),
        -311, -313, -311, 2*datum(-321);
      data.values, 2*datum(-311), 111, 221, 331, 113, 223, 211, 111;
      data.values, 211, 111, -311, -313, -311, 2*datum(-321), 2*datum(-311),
        111;
      data.values, 221, 331, 113, 223, 211, 111, 211, 111;
      data.values, 20*datum(0), 3*datum(111), 2*datum(221), 331, 113,
        223, 3*datum(211), -211;
      data.values, 111, -211, 111, 211, 111, 211, -211, 111;
      data.values, 113, 111, 223, 2*datum(111), -311, 4*datum(211),
        2*datum(111), 2*datum(211);
      data.values, 111, 7*datum(211), 7*datum(111), 113, 221, 2*datum(223),
        2*datum(-211), -213;
      data.values, 4*datum(-211), -213, -211, -213, -211, 2*datum(211),
        2, 2*datum(0);
      data.values, -321, -323, -311, -321, -311, 2*datum(-321), -211, -213;
      data.values, 2*datum(-211), 211, -321, -323, -311, -321, -311,
        2*datum(-321);
      data.values, -211, -213, 2*datum(-211), 211, 46*datum(0),
        3*datum(111), 113, 2*datum(221);
      data.values, 331, 2*datum(223), -311, 3*datum(-211), -213,
        8*datum(111), 113, 3*datum(211);
      data.values, 213, 2*datum(111), -211, 3*datum(111), 113, 111,
        2*datum(113), 221;
      data.values, 331, 223, 111, 221, 331, 113, 223, 113;
      data.values, 2*datum(223), 2*datum(221), 3*datum(111), 221,
        113, 223, 4*datum(211), 3*datum(-211);
      data.values, -213, -211, 5*datum(111), -321, 3*datum(211),
        3*datum(111), 2*datum(211), 2*datum(111);
      data.values, 2*datum(-211), -213, 3*datum(111), 221, 113, 223,
        6*datum(111), 3*datum(0);
      data.values, 221, 331, 333, 321, 311, 221, 331, 333;
      data.values, 321, 311, 19*datum(0), 3, 5*datum(0), -11, 0, 2*datum(111);
      data.values, -211, -11, 11, 2*datum(221), 3*datum(0), 111,
        22*datum(0), 111;
      data.values, 2*datum(0), 22, 111, 5*datum(0), 111, 12*datum(0),
        2*datum(21), 11*datum(0);
      data.values, 2*datum(21), 2*datum(-6), 111*datum(0), -211,
        2*datum(111), -211, 3*datum(111), -211;
      data.values, 111, 211, 15*datum(0), 111, 6*datum(0), 111, -211,
        9*datum(0);
      data.values, 111, -211, 9*datum(0), 111, -211, 111, -211, 4*datum(0);
      data.values, 111, -211, 111, -211, 4*datum(0), -211, 4*datum(0), 111;
      data.values, -211, 111, -211, 4*datum(0), 111, -211, 111, -211;
      data.values, 4*datum(0), -211, 3*datum(0), -211, 5*datum(0),
        111, 211, 3*datum(0);
      data.values, 111, 10*datum(0), 2*datum(111), 211, -211, 211, -211;
      FEM_DO_SAFE(i, 1, 918) {
        data, kfdp(i, 3);
      }
    }
    {
      fem::data_values data;
      data.values, 7*datum(0), 2212, 3122, 3212, 3214, 2112, 2114, 2212;
      data.values, 2112, 3122, 3212, 3214, 2112, 2114, 2212, 2112;
      data.values, 50*datum(0), 3*datum(3), 1, 12*datum(0), 2112,
        43*datum(0), 3322, 949*datum(0);
      FEM_DO_SAFE(i, 919, 2000) {
        data, kfdp(i, 3);
      }
    }
    {
      fem::data_values data;
      data.values, 83*datum(0), 3*datum(111), 9*datum(0), -211,
        3*datum(0), 111, 2*datum(-211), 0;
      data.values, 111, 0, 2*datum(111), 113, 221, 111, -213, -211;
      data.values, 211, 123*datum(0), 13*datum(81), 37*datum(0), 111,
        3*datum(211), 111, 5*datum(0);
      data.values, -211, 111, -211, 111, 2*datum(0), 111, 3*datum(211), 111;
      data.values, 5*datum(0), -211, 111, -211, 111, 50*datum(0),
        2*datum(111), 2*datum(-211);
      data.values, 2*datum(111), -211, 211, 3*datum(111), 211,
        14*datum(111), 221, 113;
      data.values, 223, 2*datum(111), 2*datum(113), 223, 2*datum(111),
        -1, 4*datum(0), -211;
      data.values, 111, -211, 211, 111, 2*datum(0), 2*datum(111),
        -211, 2*datum(0);
      data.values, -211, 111, -211, 211, 111, 2*datum(0), 2*datum(111), -211;
      data.values, 96*datum(0), 6*datum(111), 3*datum(-211), -213,
        4*datum(111), 113, 6*datum(111), 3*datum(-211);
      data.values, 3*datum(111), 2*datum(-211), 2*datum(111), 3*datum(-211),
        12*datum(111), 6*datum(0), -321, -311;
      data.values, 3*datum(0), -321, -311, 19*datum(0), -3, 11*datum(0),
        -11, 280*datum(0);
      data.values, 111, -211, 3*datum(0), 111, 29*datum(0), -211,
        111, 5*datum(0);
      data.values, -211, 111, 50*datum(0), 2101, 2103, 2*datum(2101),
        1006*datum(0);
      FEM_DO_SAFE(i, 1, 2000) {
        data, kfdp(i, 4);
      }
    }
    {
      fem::data_values data;
      data.values, 85*datum(0), 111, 15*datum(0), 111, 7*datum(0),
        111, 0, 2*datum(111);
      data.values, 175*datum(0), 111, -211, 111, 7*datum(0), 2*datum(111),
        4*datum(0), 111;
      data.values, -211, 111, 7*datum(0), 2*datum(111), 93*datum(0),
        111, -211, 111;
      data.values, 3*datum(0), 111, -211, 4*datum(0), 111, -211, 111,
        3*datum(0);
      data.values, 111, -211, 1571*datum(0);
      FEM_DO_SAFE(i, 1, 2000) {
        data, kfdp(i, 5);
      }
    }
    {
      fem::data_values data;
      data.values, "d", "u", "s", "c", "b", "t", "l", "h";
      data.values, 2*datum(" "), "e", "nu_e", "mu", "nu_mu", "tau",
        "nu_tau", "chi";
      data.values, "nu_chi", 2*datum(" "), "g", "gamma", "Z", "W",
        "H", 6*datum(" ");
      data.values, "Z'", "Z\"", "W'", "H'", "H\"", "H", 2*datum(" "), "R";
      data.values, 40*datum(" "), "specflav", "rndmflav", "phasespa",
        "c-hadron", "b-hadron", "t-hadron", "l-hadron";
      data.values, "h-hadron", "Wvirt", "diquark", "cluster",
        "string", "indep.", "CMshower", "SPHEaxis";
      data.values, "THRUaxis", "CLUSjet", "CELLjet", "table", " ",
        "pi", 2*datum("K"), 2*datum("D");
      data.values, "D_s", 2*datum("B"), "B_s", " ", "pi", "eta",
        "eta'", "eta_c";
      data.values, "eta_b", "eta_t", "eta_l", "eta_h", 2*datum(" "),
        "rho", 2*datum("K*"), 2*datum("D*");
      data.values, "D*_s", 2*datum("B*"), "B*_s", " ", "rho",
        "omega", "phi", "J/psi";
      data.values, "Upsilon", "Theta", "Theta_l", "Theta_h", 2*datum(" "),
        "b_1", 2*datum("K_1"), 2*datum("D_1");
      data.values, "D_1s", 2*datum("B_1"), "B_1s", " ", "b_1", "h_1",
        "h'_1", "h_1c";
      data.values, "h_1b", "h_1t", "h_1l", "h_1h", 2*datum(" "),
        "a_0", 2*datum("K*_0"), 2*datum("D*_0");
      data.values, "D*_0s", 2*datum("B*_0"), "B*_0s", " ", "a_0",
        "f_0", "f'_0", "chi_0c";
      data.values, "chi_0b", "chi_0t", "chi_0l", "chi_0h", 2*datum(" "),
        "a_1", 2*datum("K*_1"), 2*datum("D*_1");
      data.values, "D*_1s", 2*datum("B*_1"), "B*_1s", " ", "a_1",
        "f_1", "f'_1", "chi_1c";
      data.values, "chi_1b", "chi_1t", "chi_1l", "chi_1h", 2*datum(" "),
        "a_2", 2*datum("K*_2"), 2*datum("D*_2");
      data.values, "D*_2s", 2*datum("B*_2"), "B*_2s", " ", "a_2",
        "f_2", "f'_2", "chi_2c";
      data.values, "chi_2b", "chi_2t", "chi_2l", "chi_2h", 2*datum(" "),
        "K_L", "K_S", 58*datum(" ");
      data.values, "pi_diffr", "n_diffr", "p_diffr", 22*datum(" "),
        "Lambda", 5*datum(" "), "Lambda_c", " ";
      data.values, 2*datum("Xi_c"), 6*datum(" "), "Lambda_b", " ",
        2*datum("Xi_b"), 6*datum(" ");
      FEM_DO_SAFE(i, 1, 331) {
        data, chaf(i);
      }
    }
    {
      fem::data_values data;
      data.values, "n", "p", " ", 3*datum("Sigma"), 2*datum("Xi"),
        " ", 3*datum("Sigma_c"), 2*datum("Xi'_c");
      data.values, "Omega_c", 4*datum(" "), 3*datum("Sigma_b"),
        2*datum("Xi'_b"), "Omega_b", 4*datum(" "), 4*datum("Delta"),
        3*datum("Sigma*");
      data.values, 2*datum("Xi*"), "Omega", 3*datum("Sigma*_c"),
        2*datum("Xi*_c"), "Omega*_c", 4*datum(" "), 3*datum("Sigma*_b"),
        2*datum("Xi*_b");
      data.values, "Omega*_b", 114*datum(" ");
      FEM_DO_SAFE(i, 332, 500) {
        data, chaf(i);
      }
    }
    {
      static const int values[] = {
        19780503, 0, 0, 97, 33, 0
      };
      fem::data_of_type<int>(FEM_VALUES_AND_SIZE),
        mrlu;
    }
  }
  //C
  //C...Purpose: to give default values to parameters and particle and
  //C...decay data.
  //C
  //C...LUDAT1, containing status codes and most parameters.
  //C
  //C...LUDAT2, with particle data and flavour treatment parameters.
  //C
  //C...LUDAT3, with particle decay parameters and data.
  //C
  //C...LUDAT4, with character strings.
  //C
  //C...LUDATR, with initial values for the random number generator.
  //C
}

struct pyinki_save
{
  arr<fem::str<26> > chalp;
  arr<fem::str<8> > chcde;
  arr<int> kcde;

  pyinki_save() :
    chalp(dimension(2), fem::fill0),
    chcde(dimension(18), fem::fill0),
    kcde(dimension(18), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
pyinki(
  common& cmn,
  str_cref chfram,
  str_cref chbeam,
  str_cref chtarg,
  float const& win)
{
  FEM_CMN_SVE(pyinki);
  common_write write(cmn);
  // COMMON lujets
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  // COMMON ludat1
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  // COMMON pypars
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  // COMMON pyint1
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  //
  // SAVE
  str_arr_ref<1> chalp(sve.chalp, dimension(2));
  str_arr_ref<1> chcde(sve.chcde, dimension(18));
  arr_ref<int> kcde(sve.kcde, dimension(18));
  //
  if (is_called_first_time) {
    {
      static const char* values[] = {
        "abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chalp;
    }
    {
      static const char* values[] = {
        "e-      ", "e+      ", "nue     ", "nue~    ", "mu-     ",
          "mu+     ", "numu    ", "numu~   ", "tau-    ", "tau+    ",
          "nutau   ", "nutau~  ", "pi+     ", "pi-     ", "n       ",
          "n~      ", "p       ", "p~      "
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chcde;
    }
    {
      static const int values[] = {
        11, -11, 12, -12, 13, -13, 14, -14, 15, -15, 16, -16, 211,
          -211, 2112, -2112, 2212, -2212
      };
      fem::data_of_type<int>(FEM_VALUES_AND_SIZE),
        kcde;
    }
  }
  //C
  //C...Identifies the two incoming particles and sets up kinematics,
  //C...including rotations and boosts to/from CM frame.
  //C
  //C...Convert character variables to lowercase and find their length.
  arr_1d<3, fem::str<8> > chcom(fem::fill0);
  chcom(1) = chfram;
  chcom(2) = chbeam;
  chcom(3) = chtarg;
  int i = fem::int0;
  arr_1d<3, int> len(fem::fill0);
  int ll = fem::int0;
  int la = fem::int0;
  arr_1d<3, fem::str<8> > chidnt(fem::fill0);
  fem::str<8> chtemp = fem::char0;
  FEM_DO_SAFE(i, 1, 3) {
    len(i) = 8;
    FEM_DOSTEP(ll, 8, 1, -1) {
      if (len(i) == ll && chcom(i)(ll, ll) == " ") {
        len(i) = ll - 1;
      }
      FEM_DO_SAFE(la, 1, 26) {
        if (chcom(i)(ll, ll) == chalp(2)(la, la)) {
          chcom(i)(ll, ll) = chalp(1)(la, la);
        }
      }
    }
    chidnt(i) = chcom(i);
    FEM_DO_SAFE(ll, 1, 6) {
      if (chidnt(i)(ll, ll + 2) == "bar") {
        chtemp = chidnt(i);
        chidnt(i) = chtemp(1, ll - 1) + str_cref("~") + chtemp(ll + 3,
          8) + str_cref("  ");
      }
    }
    FEM_DO_SAFE(ll, 1, 8) {
      if (chidnt(i)(ll, ll) == "_") {
        chtemp = chidnt(i);
        chidnt(i) = chtemp(1, ll - 1) + chtemp(ll + 1, 8) + str_cref(" ");
      }
    }
  }
  //C
  //C...Set initial state. Error for unknown codes. Reset variables.
  cmn.n = 2;
  int j = fem::int0;
  FEM_DO_SAFE(i, 1, 2) {
    k(i, 2) = 0;
    FEM_DO_SAFE(j, 1, 18) {
      if (chidnt(i + 1) == chcde(j)) {
        k(i, 2) = kcde(j);
      }
    }
    p(i, 5) = ulmass(cmn, k(i, 2));
    mint(40 + i) = 1;
    if (fem::iabs(k(i, 2)) > 100) {
      mint(40 + i) = 2;
    }
    FEM_DO_SAFE(j, 1, 5) {
      v(i, j) = 0.f;
    }
  }
  if (k(1, 2) == 0) {
    write(mstu(11),
      "(1x,'Error: unrecognized beam particle ''',a,'''.',/,1x,"
      "'Execution stopped!')"),
      chbeam(1, len(2));
  }
  if (k(2, 2) == 0) {
    write(mstu(11),
      "(1x,'Error: unrecognized target particle ''',a,'''.',/,1x,"
      "'Execution stopped!')"),
      chtarg(1, len(3));
  }
  if (k(1, 2) == 0 || k(2, 2) == 0) {
    FEM_STOP(0);
  }
  FEM_DO_SAFE(j, 6, 10) {
    vint(j) = 0.f;
  }
  fem::str<76> chinit = " ";
  //C
  //C...Set up kinematics for events defined in CM frame.
  int loffs = fem::int0;
  float s = fem::float0;
  if (chcom(1)(1, 2) == "cm") {
    if (chcom(2)(1, 1) != "e") {
      loffs = (34 - (len(2) + len(3))) / 2;
      chinit(loffs + 1, 76) = "PYTHIA will be initialized for a " + chcom(2)(1,
        len(2)) + str_cref("-") + chcom(3)(1, len(3)) + str_cref(
        " collider") + str_cref(" ");
    }
    else {
      loffs = (33 - (len(2) + len(3))) / 2;
      chinit(loffs + 1, 76) = "PYTHIA will be initialized for an " +
        chcom(2)(1, len(2)) + str_cref("-") + chcom(3)(1, len(3)) +
        str_cref(" collider") + str_cref(" ");
    }
    //C        WRITE(MSTU(11),1200) CHINIT
    //C        WRITE(MSTU(11),1300) WIN
    s = fem::pow2(win);
    p(1, 1) = 0.f;
    p(1, 2) = 0.f;
    p(2, 1) = 0.f;
    p(2, 2) = 0.f;
    p(1, 3) = fem::sqrt((fem::pow2((s - fem::pow2(p(1, 5)) - fem::pow2(p(2,
      5)))) - fem::pow2((2.f * p(1, 5) * p(2, 5)))) / (4.f * s));
    p(2, 3) = -p(1, 3);
    p(1, 4) = fem::sqrt(fem::pow2(p(1, 3)) + fem::pow2(p(1, 5)));
    p(2, 4) = fem::sqrt(fem::pow2(p(2, 3)) + fem::pow2(p(2, 5)));
    //C
    //C...Set up kinematics for fixed target events.
  }
  else if (chcom(1)(1, 3) == "fix") {
    loffs = (29 - (len(2) + len(3))) / 2;
    chinit(loffs + 1, 76) = "PYTHIA will be initialized for " + chcom(2)(1,
      len(2)) + str_cref(" on ") + chcom(3)(1, len(3)) + str_cref(
      " fixed target") + str_cref(" ");
    //C        WRITE(MSTU(11),1200) CHINIT
    //C        WRITE(MSTU(11),1400) WIN
    p(1, 1) = 0.f;
    p(1, 2) = 0.f;
    p(2, 1) = 0.f;
    p(2, 2) = 0.f;
    p(1, 3) = win;
    p(1, 4) = fem::sqrt(fem::pow2(p(1, 3)) + fem::pow2(p(1, 5)));
    p(2, 3) = 0.f;
    p(2, 4) = p(2, 5);
    s = fem::pow2(p(1, 5)) + fem::pow2(p(2, 5)) + 2.f * p(2, 4) * p(1, 4);
    vint(10) = p(1, 3) / (p(1, 4) + p(2, 4));
    lurobo(cmn, 0.f, 0.f, 0.f, 0.f, -vint(10));
    //C        WRITE(MSTU(11),1500) SQRT(S)
    //C
    //C...Set up kinematics for events in user-defined frame.
  }
  else if (chcom(1)(1, 3) == "use") {
    loffs = (13 - (len(1) + len(2))) / 2;
    chinit(loffs + 1, 76) = "PYTHIA will be initialized for " + chcom(2)(1,
      len(2)) + str_cref(" on ") + chcom(3)(1, len(3)) + str_cref(
      "user-specified configuration") + str_cref(" ");
    //C        WRITE(MSTU(11),1200) CHINIT
    //C        WRITE(MSTU(11),1600)
    //C        WRITE(MSTU(11),1700) CHCOM(2),P(1,1),P(1,2),P(1,3)
    //C        WRITE(MSTU(11),1700) CHCOM(3),P(2,1),P(2,2),P(2,3)
    p(1, 4) = fem::sqrt(fem::pow2(p(1, 1)) + fem::pow2(p(1, 2)) +
      fem::pow2(p(1, 3)) + fem::pow2(p(1, 5)));
    p(2, 4) = fem::sqrt(fem::pow2(p(2, 1)) + fem::pow2(p(2, 2)) +
      fem::pow2(p(2, 3)) + fem::pow2(p(2, 5)));
    FEM_DO_SAFE(j, 1, 3) {
      vint(7 + j) = fem::sngl((fem::dble(p(1, j)) + fem::dble(p(2,
        j))) / fem::dble(p(1, 4) + p(2, 4)));
    }
    lurobo(cmn, 0.f, 0.f, -vint(8), -vint(9), -vint(10));
    vint(7) = ulangl(cmn, p(1, 1), p(1, 2));
    lurobo(cmn, 0.f, -vint(7), 0.f, 0.f, 0.f);
    vint(6) = ulangl(cmn, p(1, 3), p(1, 1));
    lurobo(cmn, -vint(6), 0.f, 0.f, 0.f, 0.f);
    s = fem::pow2(p(1, 5)) + fem::pow2(p(2, 5)) + 2.f * (p(1, 4) * p(2,
      4) - p(1, 3) * p(2, 3));
    //C        WRITE(MSTU(11),1500) SQRT(S)
    //C
    //C...Unknown frame. Error for too low CM energy.
  }
  else {
    write(mstu(11),
      "(1x,'Error: unrecognized coordinate frame ''',a,'''.',/,1x,"
      "'Execution stopped!')"),
      chfram(1, len(1));
    FEM_STOP(0);
  }
  if (s < fem::pow2(parp(2))) {
    write(mstu(11),
      "(1x,'Error: too low CM energy,',f8.3,' GeV for event ','generation.',/,"
      "1x,'Execution stopped!')"),
      fem::sqrt(s);
    FEM_STOP(0);
  }
  //C
  //C...Save information on incoming particles.
  mint(11) = k(1, 2);
  mint(12) = k(2, 2);
  mint(43) = 2 * mint(41) + mint(42) - 2;
  vint(1) = fem::sqrt(s);
  vint(2) = s;
  vint(3) = p(1, 5);
  vint(4) = p(2, 5);
  vint(5) = p(1, 3);
  //C
  //C...Store constants to be used in generation.
  if (mstp(82) <= 1) {
    vint(149) = 4.f * fem::pow2(parp(81)) / s;
  }
  if (mstp(82) >= 2) {
    vint(149) = 4.f * fem::pow2(parp(82)) / s;
  }
  //C
  //C...Formats for initialization and error information.
  //Clin 1200 FORMAT(/1X,78('=')/1X,'I',76X,'I'/1X,'I',A76,'I')
  //C 1300 FORMAT(1X,'I',18X,'at',1X,F10.3,1X,'GeV center-of-mass energy',
  //C     &19X,'I'/1X,'I',76X,'I'/1X,78('='))
  //C 1400 FORMAT(1X,'I',22X,'at',1X,F10.3,1X,'GeV/c lab-momentum',22X,'I')
  //C 1500 FORMAT(1X,'I',76X,'I'/1X,'I',11X,'corresponding to',1X,F10.3,1X,
  //C     &'GeV center-of-mass energy',12X,'I'/1X,'I',76X,'I'/1X,78('='))
  //C 1600 FORMAT(1X,'I',76X,'I'/1X,'I',24X,'px (GeV/c)',3X,'py (GeV/c)',3X,
  //C     &'pz (GeV/c)',16X,'I')
  //Clin 1700 FORMAT(1X,'I',15X,A8,3(2X,F10.3,1X),15X,'I')
  //C
}

//C
//C*********************************************************************
//C
void
pywidt(
  common& cmn,
  int const& kflr,
  float const& rmas,
  arr_ref<float> wdtp,
  arr_ref<float, 2> wdte)
{
  wdtp(dim1(0, 40));
  wdte(dim1(0, 40).dim2(0, 5));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<float, 2> vckm(cmn.vckm, dimension(4, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_cref<int, 2> kfdp(cmn.kfdp, dimension(2000, 5));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<float, 2> wids(cmn.wids, dim1(21, 40).dim2(3));
  //
  int kfla = fem::int0;
  float sqm = fem::float0;
  float as = fem::float0;
  float aem = fem::float0;
  float xw = fem::float0;
  float radc = fem::float0;
  int i = fem::int0;
  int j = fem::int0;
  int idc = fem::int0;
  float rm1 = fem::float0;
  float rm2 = fem::float0;
  float wid2 = fem::float0;
  float ei = fem::float0;
  float ai = fem::float0;
  float vi = fem::float0;
  float sqmz = fem::float0;
  float gzmz = fem::float0;
  float ggi = fem::float0;
  float gzi = fem::float0;
  float zzi = fem::float0;
  float ef = fem::float0;
  float af = fem::float0;
  float vf = fem::float0;
  float ggf = fem::float0;
  float gzf = fem::float0;
  float zzf = fem::float0;
  float cf = fem::float0;
  float etare = fem::float0;
  float etaim = fem::float0;
  float eps = fem::float0;
  float root = fem::float0;
  float rln = fem::float0;
  float phire = fem::float0;
  float phiim = fem::float0;
  float eta2 = fem::float0;
  float ej = fem::float0;
  int jl = fem::int0;
  float aj = fem::float0;
  float vj = fem::float0;
  float epsp = fem::float0;
  float psire = fem::float0;
  float psiim = fem::float0;
  float phirep = fem::float0;
  float phiimp = fem::float0;
  float psirep = fem::float0;
  float psiimp = fem::float0;
  float fxyre = fem::float0;
  float fxyim = fem::float0;
  float f1re = fem::float0;
  float f1im = fem::float0;
  float api = fem::float0;
  float vpi = fem::float0;
  float sqmzp = fem::float0;
  float gzpmzp = fem::float0;
  float gzpi = fem::float0;
  float zzpi = fem::float0;
  float zpzpi = fem::float0;
  float apf = fem::float0;
  float vpf = fem::float0;
  float gzpf = fem::float0;
  float zzpf = fem::float0;
  float zpzpf = fem::float0;
  //C
  //C...Calculates full and partial widths of resonances.
  //C
  //C...Some common constants.
  kfla = fem::iabs(kflr);
  sqm = fem::pow2(rmas);
  as = ulalps(cmn, sqm);
  aem = paru(101);
  xw = paru(102);
  radc = 1.f + as / paru(1);
  //C
  //C...Reset width information.
  FEM_DO_SAFE(i, 0, 40) {
    wdtp(i) = 0.f;
    FEM_DO_SAFE(j, 0, 5) {
      wdte(i, j) = 0.f;
    }
  }
  //C
  if (kfla == 21) {
    //C...QCD:
    FEM_DO_SAFE(i, 1, mdcy(21, 3)) {
      idc = i + mdcy(21, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_110;
      }
      if (i <= 8) {
        //C...QCD -> q + qb
        wdtp(i) = (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f, 1.f - 4.f * rm1));
        wid2 = 1.f;
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
      }
      statement_110:;
    }
    //C
  }
  else if (kfla == 23) {
    //C...Z0:
    if (mint(61) == 1) {
      ei = kchg(fem::iabs(mint(15)), 1) / 3.f;
      ai = fem::sign(1.f, ei);
      vi = ai - 4.f * ei * xw;
      sqmz = fem::pow2(pmas(23, 1));
      gzmz = pmas(23, 2) * pmas(23, 1);
      ggi = fem::pow2(ei);
      gzi = ei * vi / (8.f * xw * (1.f - xw)) * sqm * (sqm - sqmz) / (
        fem::pow2((sqm - sqmz)) + fem::pow2(gzmz));
      zzi = (fem::pow2(vi) + fem::pow2(ai)) / fem::pow2((16.f * xw * (
        1.f - xw))) * fem::pow2(sqm) / (fem::pow2((sqm - sqmz)) +
        fem::pow2(gzmz));
      if (mstp(43) == 1) {
        //C...Only gamma* production included
        gzi = 0.f;
        zzi = 0.f;
      }
      else if (mstp(43) == 2) {
        //C...Only Z0 production included
        ggi = 0.f;
        gzi = 0.f;
      }
    }
    else if (mint(61) == 2) {
      vint(111) = 0.f;
      vint(112) = 0.f;
      vint(114) = 0.f;
    }
    FEM_DO_SAFE(i, 1, mdcy(23, 3)) {
      idc = i + mdcy(23, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_120;
      }
      if (i <= 8) {
        //C...Z0 -> q + qb
        ef = kchg(i, 1) / 3.f;
        af = fem::sign(1.f, ef + 0.1f);
        vf = af - 4.f * ef * xw;
        if (mint(61) == 0) {
          wdtp(i) = 3.f * (fem::pow2(vf) * (1.f + 2.f * rm1) +
            fem::pow2(af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
        }
        else if (mint(61) == 1) {
          wdtp(i) = 3.f * ((ggi * fem::pow2(ef) + gzi * ef * vf +
            zzi * fem::pow2(vf)) * (1.f + 2.f * rm1) + zzi *
            fem::pow2(af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
        }
        else if (mint(61) == 2) {
          ggf = 3.f * fem::pow2(ef) * (1.f + 2.f * rm1) * fem::sqrt(
            fem::max(0.f, 1.f - 4.f * rm1)) * radc;
          gzf = 3.f * ef * vf * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
          zzf = 3.f * (fem::pow2(vf) * (1.f + 2.f * rm1) + fem::pow2(
            af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
            4.f * rm1)) * radc;
        }
        wid2 = 1.f;
      }
      else if (i <= 16) {
        //C...Z0 -> l+ + l-, nu + nub
        ef = kchg(i + 2, 1) / 3.f;
        af = fem::sign(1.f, ef + 0.1f);
        vf = af - 4.f * ef * xw;
        wdtp(i) = (fem::pow2(vf) * (1.f + 2.f * rm1) + fem::pow2(
          af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
          4.f * rm1));
        if (mint(61) == 0) {
          wdtp(i) = (fem::pow2(vf) * (1.f + 2.f * rm1) + fem::pow2(
            af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
            4.f * rm1));
        }
        else if (mint(61) == 1) {
          wdtp(i) = ((ggi * fem::pow2(ef) + gzi * ef * vf + zzi *
            fem::pow2(vf)) * (1.f + 2.f * rm1) + zzi * fem::pow2(
            af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
            4.f * rm1));
        }
        else if (mint(61) == 2) {
          ggf = fem::pow2(ef) * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          gzf = ef * vf * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          zzf = (fem::pow2(vf) * (1.f + 2.f * rm1) + fem::pow2(af) * (
            1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f - 4.f *
            rm1));
        }
        wid2 = 1.f;
      }
      else {
        //C...Z0 -> H+ + H-
        cf = 2.f * (1.f - 2.f * xw);
        if (mint(61) == 0) {
          wdtp(i) = 0.25f * fem::pow2(cf) * (1.f - 4.f * rm1) *
            fem::sqrt(fem::max(0.f, 1.f - 4.f * rm1));
        }
        else if (mint(61) == 1) {
          wdtp(i) = 0.25f * (ggi + gzi * cf + zzi * fem::pow2(cf)) * (
            1.f - 4.f * rm1) * fem::sqrt(fem::max(0.f, 1.f - 4.f *
            rm1));
        }
        else if (mint(61) == 2) {
          ggf = 0.25f * (1.f - 4.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          gzf = 0.25f * cf * (1.f - 4.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          zzf = 0.25f * fem::pow2(cf) * (1.f - 4.f * rm1) * fem::sqrt(
            fem::max(0.f, 1.f - 4.f * rm1));
        }
        wid2 = wids(37, 1);
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
        //Clin-4/2008 modified a la pythia6115.f to avoid undefined values (GGF,GZF,ZZF):
        //C          VINT(111)=VINT(111)+GGF*WID2
        //C          VINT(112)=VINT(112)+GZF*WID2
        //C          VINT(114)=VINT(114)+ZZF*WID2
        if (mint(61) == 2) {
          vint(111) += ggf * wid2;
          vint(112) += gzf * wid2;
          vint(114) += zzf * wid2;
        }
        //Clin-4/2008-end
      }
      statement_120:;
    }
    if (mstp(43) == 1) {
      //C...Only gamma* production included
      vint(112) = 0.f;
      vint(114) = 0.f;
    }
    else if (mstp(43) == 2) {
      //C...Only Z0 production included
      vint(111) = 0.f;
      vint(112) = 0.f;
    }
    //C
  }
  else if (kfla == 24) {
    //C...W+/-:
    FEM_DO_SAFE(i, 1, mdcy(24, 3)) {
      idc = i + mdcy(24, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_130;
      }
      if (i <= 16) {
        //C...W+/- -> q + qb'
        wdtp(i) = 3.f * (2.f - rm1 - rm2 - fem::pow2((rm1 - rm2))) *
          fem::sqrt(fem::max(0.f, fem::pow2((1.f - rm1 - rm2)) -
          4.f * rm1 * rm2)) * vckm((i - 1) / 4 + 1, fem::mod(i - 1,
          4) + 1) * radc;
        wid2 = 1.f;
      }
      else {
        //C...W+/- -> l+/- + nu
        wdtp(i) = (2.f - rm1 - rm2 - fem::pow2((rm1 - rm2))) *
          fem::sqrt(fem::max(0.f, fem::pow2((1.f - rm1 - rm2)) -
          4.f * rm1 * rm2));
        wid2 = 1.f;
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
      }
      statement_130:;
    }
    //C
  }
  else if (kfla == 25) {
    //C...H0:
    FEM_DO_SAFE(i, 1, mdcy(25, 3)) {
      idc = i + mdcy(25, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_170;
      }
      if (i <= 8) {
        //C...H0 -> q + qb
        wdtp(i) = 3.f * rm1 * (1.f - 4.f * rm1) * fem::sqrt(fem::max(0.f,
          1.f - 4.f * rm1)) * radc;
        wid2 = 1.f;
      }
      else if (i <= 12) {
        //C...H0 -> l+ + l-
        wdtp(i) = rm1 * (1.f - 4.f * rm1) * fem::sqrt(fem::max(0.f,
          1.f - 4.f * rm1));
        wid2 = 1.f;
      }
      else if (i == 13) {
        //C...H0 -> g + g; quark loop contribution only
        etare = 0.f;
        etaim = 0.f;
        FEM_DO_SAFE(j, 1, 2 * mstp(1)) {
          eps = fem::pow2((2.f * pmas(j, 1) / rmas));
          if (eps <= 1.f) {
            if (eps > 1.e-4f) {
              root = fem::sqrt(1.f - eps);
              rln = fem::log((1.f + root) / (1.f - root));
            }
            else {
              rln = fem::log(4.f / eps - 2.f);
            }
            phire = 0.25f * (fem::pow2(rln) - fem::pow2(paru(1)));
            phiim = 0.5f * paru(1) * rln;
          }
          else {
            phire = -fem::pow2((fem::asin(1.f / fem::sqrt(eps))));
            phiim = 0.f;
          }
          etare += 0.5f * eps * (1.f + (eps - 1.f) * phire);
          etaim += 0.5f * eps * (eps - 1.f) * phiim;
        }
        eta2 = fem::pow2(etare) + fem::pow2(etaim);
        wdtp(i) = fem::pow2((as / paru(1))) * eta2;
        wid2 = 1.f;
      }
      else if (i == 14) {
        //C...H0 -> gamma + gamma; quark, charged lepton and W loop contributions
        etare = 0.f;
        etaim = 0.f;
        FEM_DO_SAFE(j, 1, 3 * mstp(1) + 1) {
          if (j <= 2 * mstp(1)) {
            ej = kchg(j, 1) / 3.f;
            eps = fem::pow2((2.f * pmas(j, 1) / rmas));
          }
          else if (j <= 3 * mstp(1)) {
            jl = 2 * (j - 2 * mstp(1)) - 1;
            ej = kchg(10 + jl, 1) / 3.f;
            eps = fem::pow2((2.f * pmas(10 + jl, 1) / rmas));
          }
          else {
            eps = fem::pow2((2.f * pmas(24, 1) / rmas));
          }
          if (eps <= 1.f) {
            if (eps > 1.e-4f) {
              root = fem::sqrt(1.f - eps);
              rln = fem::log((1.f + root) / (1.f - root));
            }
            else {
              rln = fem::log(4.f / eps - 2.f);
            }
            phire = 0.25f * (fem::pow2(rln) - fem::pow2(paru(1)));
            phiim = 0.5f * paru(1) * rln;
          }
          else {
            phire = -fem::pow2((fem::asin(1.f / fem::sqrt(eps))));
            phiim = 0.f;
          }
          if (j <= 2 * mstp(1)) {
            etare += 0.5f * 3.f * fem::pow2(ej) * eps * (1.f + (eps -
              1.f) * phire);
            etaim += 0.5f * 3.f * fem::pow2(ej) * eps * (eps - 1.f) * phiim;
          }
          else if (j <= 3 * mstp(1)) {
            etare += 0.5f * fem::pow2(ej) * eps * (1.f + (eps - 1.f) * phire);
            etaim += 0.5f * fem::pow2(ej) * eps * (eps - 1.f) * phiim;
          }
          else {
            etare = etare - 0.5f - 0.75f * eps * (1.f + (eps - 2.f) * phire);
            etaim += 0.75f * eps * (eps - 2.f) * phiim;
          }
        }
        eta2 = fem::pow2(etare) + fem::pow2(etaim);
        wdtp(i) = fem::pow2((aem / paru(1))) * 0.5f * eta2;
        wid2 = 1.f;
      }
      else if (i == 15) {
        //C...H0 -> gamma + Z0; quark, charged lepton and W loop contributions
        etare = 0.f;
        etaim = 0.f;
        FEM_DO_SAFE(j, 1, 3 * mstp(1) + 1) {
          if (j <= 2 * mstp(1)) {
            ej = kchg(j, 1) / 3.f;
            aj = fem::sign(1.f, ej + 0.1f);
            vj = aj - 4.f * ej * xw;
            eps = fem::pow2((2.f * pmas(j, 1) / rmas));
            epsp = fem::pow2((2.f * pmas(j, 1) / pmas(23, 1)));
          }
          else if (j <= 3 * mstp(1)) {
            jl = 2 * (j - 2 * mstp(1)) - 1;
            ej = kchg(10 + jl, 1) / 3.f;
            aj = fem::sign(1.f, ej + 0.1f);
            vj = ai - 4.f * ej * xw;
            eps = fem::pow2((2.f * pmas(10 + jl, 1) / rmas));
            epsp = fem::pow2((2.f * pmas(10 + jl, 1) / pmas(23, 1)));
          }
          else {
            eps = fem::pow2((2.f * pmas(24, 1) / rmas));
            epsp = fem::pow2((2.f * pmas(24, 1) / pmas(23, 1)));
          }
          if (eps <= 1.f) {
            root = fem::sqrt(1.f - eps);
            if (eps > 1.e-4f) {
              rln = fem::log((1.f + root) / (1.f - root));
            }
            else {
              rln = fem::log(4.f / eps - 2.f);
            }
            phire = 0.25f * (fem::pow2(rln) - fem::pow2(paru(1)));
            phiim = 0.5f * paru(1) * rln;
            psire = -(1.f + 0.5f * root * rln);
            psiim = 0.5f * paru(1) * root;
          }
          else {
            phire = -fem::pow2((fem::asin(1.f / fem::sqrt(eps))));
            phiim = 0.f;
            psire = -(1.f + fem::sqrt(eps - 1.f) * fem::asin(1.f /
              fem::sqrt(eps)));
            psiim = 0.f;
          }
          if (epsp <= 1.f) {
            root = fem::sqrt(1.f - epsp);
            if (epsp > 1.e-4f) {
              rln = fem::log((1.f + root) / (1.f - root));
            }
            else {
              rln = fem::log(4.f / epsp - 2.f);
            }
            phirep = 0.25f * (fem::pow2(rln) - fem::pow2(paru(1)));
            phiimp = 0.5f * paru(1) * rln;
            psirep = -(1.f + 0.5f * root * rln);
            psiimp = 0.5f * paru(1) * root;
          }
          else {
            phirep = -fem::pow2((fem::asin(1.f / fem::sqrt(epsp))));
            phiimp = 0.f;
            psirep = -(1.f + fem::sqrt(epsp - 1.f) * fem::asin(1.f /
              fem::sqrt(epsp)));
            psiimp = 0.f;
          }
          fxyre = eps * epsp / (8.f * (eps - epsp)) * (1.f - eps *
            epsp / (eps - epsp) * (phire - phirep) + 2.f * eps / (
            eps - epsp) * (psire - psirep));
          fxyim = eps * epsp / (8.f * (eps - epsp)) * (-eps * epsp / (
            eps - epsp) * (phiim - phiimp) + 2.f * eps / (eps -
            epsp) * (psiim - psiimp));
          f1re = eps * epsp / (2.f * (eps - epsp)) * (phire - phirep);
          f1im = eps * epsp / (2.f * (eps - epsp)) * (phiim - phiimp);
          if (j <= 2 * mstp(1)) {
            etare = etare - 3.f * ej * vj * (fxyre - 0.25f * f1re);
            etaim = etaim - 3.f * ej * vj * (fxyim - 0.25f * f1im);
          }
          else if (j <= 3 * mstp(1)) {
            etare = etare - ej * vj * (fxyre - 0.25f * f1re);
            etaim = etaim - ej * vj * (fxyim - 0.25f * f1im);
          }
          else {
            etare = etare - fem::sqrt(1.f - xw) * (((1.f + 2.f /
              eps) * xw / fem::sqrt(1.f - xw) - (5.f + 2.f / eps)) *
              fxyre + (3.f - xw / fem::sqrt(1.f - xw)) * f1re);
            etaim = etaim - fem::sqrt(1.f - xw) * (((1.f + 2.f /
              eps) * xw / fem::sqrt(1.f - xw) - (5.f + 2.f / eps)) *
              fxyim + (3.f - xw / fem::sqrt(1.f - xw)) * f1im);
          }
        }
        eta2 = fem::pow2(etare) + fem::pow2(etaim);
        wdtp(i) = fem::pow2((aem / paru(1))) * fem::pow3((1.f -
          fem::pow2((pmas(23, 1) / rmas)))) / xw * eta2;
        wid2 = wids(23, 2);
      }
      else {
        //C...H0 -> Z0 + Z0, W+ + W-
        wdtp(i) = (1.f - 4.f * rm1 + 12.f * fem::pow2(rm1)) *
          fem::sqrt(fem::max(0.f, 1.f - 4.f * rm1)) / (2.f * (18 -
          i));
        wid2 = wids(7 + i, 1);
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
      }
      statement_170:;
    }
    //C
  }
  else if (kfla == 32) {
    //C...Z'0:
    if (mint(61) == 1) {
      ei = kchg(fem::iabs(mint(15)), 1) / 3.f;
      ai = fem::sign(1.f, ei);
      vi = ai - 4.f * ei * xw;
      sqmz = fem::pow2(pmas(23, 1));
      gzmz = pmas(23, 2) * pmas(23, 1);
      api = fem::sign(1.f, ei);
      vpi = api - 4.f * ei * xw;
      sqmzp = fem::pow2(pmas(32, 1));
      gzpmzp = pmas(32, 2) * pmas(32, 1);
      ggi = fem::pow2(ei);
      gzi = ei * vi / (8.f * xw * (1.f - xw)) * sqm * (sqm - sqmz) / (
        fem::pow2((sqm - sqmz)) + fem::pow2(gzmz));
      gzpi = ei * vpi / (8.f * xw * (1.f - xw)) * sqm * (sqm -
        sqmzp) / (fem::pow2((sqm - sqmzp)) + fem::pow2(gzpmzp));
      zzi = (fem::pow2(vi) + fem::pow2(ai)) / fem::pow2((16.f * xw * (
        1.f - xw))) * fem::pow2(sqm) / (fem::pow2((sqm - sqmz)) +
        fem::pow2(gzmz));
      zzpi = 2.f * (vi * vpi + ai * api) / fem::pow2((16.f * xw * (
        1.f - xw))) * fem::pow2(sqm) * ((sqm - sqmz) * (sqm -
        sqmzp) + gzmz * gzpmzp) / ((fem::pow2((sqm - sqmz)) +
        fem::pow2(gzmz)) * (fem::pow2((sqm - sqmzp)) + fem::pow2(
        gzpmzp)));
      zpzpi = (fem::pow2(vpi) + fem::pow2(api)) / fem::pow2((16.f *
        xw * (1.f - xw))) * fem::pow2(sqm) / (fem::pow2((sqm -
        sqmzp)) + fem::pow2(gzpmzp));
      if (mstp(44) == 1) {
        //C...Only gamma* production included
        gzi = 0.f;
        gzpi = 0.f;
        zzi = 0.f;
        zzpi = 0.f;
        zpzpi = 0.f;
      }
      else if (mstp(44) == 2) {
        //C...Only Z0 production included
        ggi = 0.f;
        gzi = 0.f;
        gzpi = 0.f;
        zzpi = 0.f;
        zpzpi = 0.f;
      }
      else if (mstp(44) == 3) {
        //C...Only Z'0 production included
        ggi = 0.f;
        gzi = 0.f;
        gzpi = 0.f;
        zzi = 0.f;
        zzpi = 0.f;
      }
      else if (mstp(44) == 4) {
        //C...Only gamma*/Z0 production included
        gzpi = 0.f;
        zzpi = 0.f;
        zpzpi = 0.f;
      }
      else if (mstp(44) == 5) {
        //C...Only gamma*/Z'0 production included
        gzi = 0.f;
        zzi = 0.f;
        zzpi = 0.f;
      }
      else if (mstp(44) == 6) {
        //C...Only Z0/Z'0 production included
        ggi = 0.f;
        gzi = 0.f;
        gzpi = 0.f;
      }
    }
    else if (mint(61) == 2) {
      vint(111) = 0.f;
      vint(112) = 0.f;
      vint(113) = 0.f;
      vint(114) = 0.f;
      vint(115) = 0.f;
      vint(116) = 0.f;
    }
    FEM_DO_SAFE(i, 1, mdcy(32, 3)) {
      idc = i + mdcy(32, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_180;
      }
      if (i <= 8) {
        //C...Z'0 -> q + qb
        ef = kchg(i, 1) / 3.f;
        af = fem::sign(1.f, ef + 0.1f);
        vf = af - 4.f * ef * xw;
        apf = fem::sign(1.f, ef + 0.1f);
        vpf = apf - 4.f * ef * xw;
        if (mint(61) == 0) {
          wdtp(i) = 3.f * (fem::pow2(vpf) * (1.f + 2.f * rm1) +
            fem::pow2(apf) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
        }
        else if (mint(61) == 1) {
          wdtp(i) = 3.f * ((ggi * fem::pow2(ef) + gzi * ef * vf +
            gzpi * ef * vpf + zzi * fem::pow2(vf) + zzpi * vf * vpf +
            zpzpi * fem::pow2(vpf)) * (1.f + 2.f * rm1) + (zzi *
            fem::pow2(af) + zzpi * af * apf + zpzpi * fem::pow2(
            apf)) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
        }
        else if (mint(61) == 2) {
          ggf = 3.f * fem::pow2(ef) * (1.f + 2.f * rm1) * fem::sqrt(
            fem::max(0.f, 1.f - 4.f * rm1)) * radc;
          gzf = 3.f * ef * vf * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
          gzpf = 3.f * ef * vpf * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
          zzf = 3.f * (fem::pow2(vf) * (1.f + 2.f * rm1) + fem::pow2(
            af) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
            4.f * rm1)) * radc;
          zzpf = 3.f * (vf * vpf * (1.f + 2.f * rm1) + af * apf * (
            1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f - 4.f *
            rm1)) * radc;
          zpzpf = 3.f * (fem::pow2(vpf) * (1.f + 2.f * rm1) +
            fem::pow2(apf) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1)) * radc;
        }
        wid2 = 1.f;
      }
      else {
        //C...Z'0 -> l+ + l-, nu + nub
        ef = kchg(i + 2, 1) / 3.f;
        af = fem::sign(1.f, ef + 0.1f);
        vf = af - 4.f * ef * xw;
        //Clin-4/2008 modified above a la pythia6115.f to avoid undefined variable API:
        //C          APF=SIGN(1.,EF+0.1)
        //C          VPF=API-4.*EF*XW
        if (i <= 10) {
          vpf = paru(127 - 2 * fem::mod(i, 2));
          apf = paru(128 - 2 * fem::mod(i, 2));
        }
        else if (i <= 12) {
          vpf = parj(186 - 2 * fem::mod(i, 2));
          apf = parj(187 - 2 * fem::mod(i, 2));
        }
        else {
          vpf = parj(194 - 2 * fem::mod(i, 2));
          apf = parj(195 - 2 * fem::mod(i, 2));
        }
        //Clin-4/2008-end
        if (mint(61) == 0) {
          wdtp(i) = (fem::pow2(vpf) * (1.f + 2.f * rm1) + fem::pow2(
            apf) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
            4.f * rm1));
        }
        else if (mint(61) == 1) {
          wdtp(i) = ((ggi * fem::pow2(ef) + gzi * ef * vf + gzpi *
            ef * vpf + zzi * fem::pow2(vf) + zzpi * vf * vpf +
            zpzpi * fem::pow2(vpf)) * (1.f + 2.f * rm1) + (zzi *
            fem::pow2(af) + zzpi * af * apf + zpzpi * fem::pow2(
            apf)) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
        }
        else if (mint(61) == 2) {
          ggf = fem::pow2(ef) * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          gzf = ef * vf * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          gzpf = ef * vpf * (1.f + 2.f * rm1) * fem::sqrt(fem::max(0.f,
            1.f - 4.f * rm1));
          zzf = (fem::pow2(vf) * (1.f + 2.f * rm1) + fem::pow2(af) * (
            1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f - 4.f *
            rm1));
          zzpf = (vf * vpf * (1.f + 2.f * rm1) + af * apf * (1.f -
            4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f - 4.f * rm1));
          zpzpf = (fem::pow2(vpf) * (1.f + 2.f * rm1) + fem::pow2(
            apf) * (1.f - 4.f * rm1)) * fem::sqrt(fem::max(0.f, 1.f -
            4.f * rm1));
        }
        wid2 = 1.f;
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
        //Clin-4/2008:
        //C          VINT(111)=VINT(111)+GGF
        //C          VINT(112)=VINT(112)+GZF
        //C          VINT(113)=VINT(113)+GZPF
        //C          VINT(114)=VINT(114)+ZZF
        //C          VINT(115)=VINT(115)+ZZPF
        //C          VINT(116)=VINT(116)+ZPZPF
        if (mint(61) == 2) {
          vint(111) += ggf;
          vint(112) += gzf;
          vint(113) += gzpf;
          vint(114) += zzf;
          vint(115) += zzpf;
          vint(116) += zpzpf;
        }
        //Clin-4/2008-end
      }
      statement_180:;
    }
    if (mstp(44) == 1) {
      //C...Only gamma* production included
      vint(112) = 0.f;
      vint(113) = 0.f;
      vint(114) = 0.f;
      vint(115) = 0.f;
      vint(116) = 0.f;
    }
    else if (mstp(44) == 2) {
      //C...Only Z0 production included
      vint(111) = 0.f;
      vint(112) = 0.f;
      vint(113) = 0.f;
      vint(115) = 0.f;
      vint(116) = 0.f;
    }
    else if (mstp(44) == 3) {
      //C...Only Z'0 production included
      vint(111) = 0.f;
      vint(112) = 0.f;
      vint(113) = 0.f;
      vint(114) = 0.f;
      vint(115) = 0.f;
    }
    else if (mstp(44) == 4) {
      //C...Only gamma*/Z0 production included
      vint(113) = 0.f;
      vint(115) = 0.f;
      vint(116) = 0.f;
    }
    else if (mstp(44) == 5) {
      //C...Only gamma*/Z'0 production included
      vint(112) = 0.f;
      vint(114) = 0.f;
      vint(115) = 0.f;
    }
    else if (mstp(44) == 6) {
      //C...Only Z0/Z'0 production included
      vint(111) = 0.f;
      vint(112) = 0.f;
      vint(113) = 0.f;
    }
    //C
  }
  else if (kfla == 37) {
    //C...H+/-:
    FEM_DO_SAFE(i, 1, mdcy(37, 3)) {
      idc = i + mdcy(37, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_190;
      }
      if (i <= 4) {
        //C...H+/- -> q + qb'
        wdtp(i) = 3.f * ((rm1 * paru(121) + rm2 / paru(121)) * (1.f -
          rm1 - rm2) - 4.f * rm1 * rm2) * fem::sqrt(fem::max(0.f,
          fem::pow2((1.f - rm1 - rm2)) - 4.f * rm1 * rm2)) * radc;
        wid2 = 1.f;
      }
      else {
        //C...H+/- -> l+/- + nu
        wdtp(i) = ((rm1 * paru(121) + rm2 / paru(121)) * (1.f - rm1 -
          rm2) - 4.f * rm1 * rm2) * fem::sqrt(fem::max(0.f, fem::pow2(
          (1.f - rm1 - rm2)) - 4.f * rm1 * rm2));
        wid2 = 1.f;
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
      }
      statement_190:;
    }
    //C
  }
  else if (kfla == 40) {
    //C...R:
    FEM_DO_SAFE(i, 1, mdcy(40, 3)) {
      idc = i + mdcy(40, 2) - 1;
      rm1 = fem::pow2((pmas(fem::iabs(kfdp(idc, 1)), 1) / rmas));
      rm2 = fem::pow2((pmas(fem::iabs(kfdp(idc, 2)), 1) / rmas));
      if (fem::sqrt(rm1) + fem::sqrt(rm2) > 1.f || mdme(idc, 1) < 0) {
        goto statement_200;
      }
      if (i <= 4) {
        //C...R -> q + qb'
        wdtp(i) = 3.f * radc;
        wid2 = 1.f;
      }
      else {
        //C...R -> l+ + l'-
        wdtp(i) = 1.f;
        wid2 = 1.f;
      }
      wdtp(0) += wdtp(i);
      if (mdme(idc, 1) > 0) {
        wdte(i, mdme(idc, 1)) = wdtp(i) * wid2;
        wdte(0, mdme(idc, 1)) += wdte(i, mdme(idc, 1));
        wdte(i, 0) = wdte(i, mdme(idc, 1));
        wdte(0, 0) += wdte(i, 0);
      }
      statement_200:;
    }
    //C
  }
  mint(61) = 0;
  //C
}

//C
//C*********************************************************************
//C
void
pyinre(
  common& cmn)
{
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_ref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_ref<float> brat(cmn.brat, dimension(2000));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<int, 2> kfpr(cmn.kfpr, dimension(200, 2));
  arr_ref<float, 2> widp(cmn.widp, dim1(21, 40).dim2(0, 40));
  arr_ref<float, 2> wide(cmn.wide, dim1(21, 40).dim2(0, 40));
  arr_ref<float, 2> wids(cmn.wids, dim1(21, 40).dim2(3));
  str_arr_ref<1> proc(cmn.proc, dim1(0, 200));
  //
  float aem = fem::float0;
  float xw = fem::float0;
  int i = fem::int0;
  int j = fem::int0;
  float wmas = fem::float0;
  float wfac = fem::float0;
  arr_1d<41, float> wdtp(dim1(0, 40), fem::fill0);
  arr_2d<41, 6, float> wdte(dim1(0, 40).dim2(0, 5), fem::fill0);
  float hcmas = fem::float0;
  float hcfac = fem::float0;
  float zmas = fem::float0;
  float zfac = fem::float0;
  float hmas = fem::float0;
  float hfac = fem::float0;
  float zpmas = fem::float0;
  float zpfac = fem::float0;
  float rmas = fem::float0;
  float rfac = fem::float0;
  int kflqm = fem::int0;
  int idc = fem::int0;
  int kc = fem::int0;
  //C
  //C...Calculates full and effective widths of guage bosons, stores masses
  //C...and widths, rescales coefficients to be used for resonance
  //C...production generation.
  //C
  //C...Calculate full and effective widths of gauge bosons.
  aem = paru(101);
  xw = paru(102);
  FEM_DO_SAFE(i, 21, 40) {
    FEM_DO_SAFE(j, 0, 40) {
      widp(i, j) = 0.f;
      wide(i, j) = 0.f;
    }
  }
  //C
  //C...W+/-:
  wmas = pmas(24, 1);
  wfac = aem / (24.f * xw) * wmas;
  pywidt(cmn, 24, wmas, wdtp, wdte);
  wids(24, 1) = ((wdte(0, 1) + wdte(0, 2)) * (wdte(0, 1) + wdte(0,
    3)) + (wdte(0, 1) + wdte(0, 2) + wdte(0, 1) + wdte(0, 3)) * (wdte(0,
    4) + wdte(0, 5)) + 2.f * wdte(0, 4) * wdte(0, 5)) / fem::pow2(
    wdtp(0));
  wids(24, 2) = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) / wdtp(0);
  wids(24, 3) = (wdte(0, 1) + wdte(0, 3) + wdte(0, 4)) / wdtp(0);
  FEM_DO_SAFE(i, 0, 40) {
    widp(24, i) = wfac * wdtp(i);
    wide(24, i) = wfac * wdte(i, 0);
  }
  //C
  //C...H+/-:
  hcmas = pmas(37, 1);
  hcfac = aem / (8.f * xw) * fem::pow2((hcmas / wmas)) * hcmas;
  pywidt(cmn, 37, hcmas, wdtp, wdte);
  wids(37, 1) = ((wdte(0, 1) + wdte(0, 2)) * (wdte(0, 1) + wdte(0,
    3)) + (wdte(0, 1) + wdte(0, 2) + wdte(0, 1) + wdte(0, 3)) * (wdte(0,
    4) + wdte(0, 5)) + 2.f * wdte(0, 4) * wdte(0, 5)) / fem::pow2(
    wdtp(0));
  wids(37, 2) = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) / wdtp(0);
  wids(37, 3) = (wdte(0, 1) + wdte(0, 3) + wdte(0, 4)) / wdtp(0);
  FEM_DO_SAFE(i, 0, 40) {
    widp(37, i) = hcfac * wdtp(i);
    wide(37, i) = hcfac * wdte(i, 0);
  }
  //C
  //C...Z0:
  zmas = pmas(23, 1);
  zfac = aem / (48.f * xw * (1.f - xw)) * zmas;
  pywidt(cmn, 23, zmas, wdtp, wdte);
  wids(23, 1) = (fem::pow2((wdte(0, 1) + wdte(0, 2))) + 2.f * (wdte(0,
    1) + wdte(0, 2)) * (wdte(0, 4) + wdte(0, 5)) + 2.f * wdte(0, 4) * wdte(0,
    5)) / fem::pow2(wdtp(0));
  wids(23, 2) = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) / wdtp(0);
  wids(23, 3) = 0.f;
  FEM_DO_SAFE(i, 0, 40) {
    widp(23, i) = zfac * wdtp(i);
    wide(23, i) = zfac * wdte(i, 0);
  }
  //C
  //C...H0:
  hmas = pmas(25, 1);
  hfac = aem / (8.f * xw) * fem::pow2((hmas / wmas)) * hmas;
  pywidt(cmn, 25, hmas, wdtp, wdte);
  wids(25, 1) = (fem::pow2((wdte(0, 1) + wdte(0, 2))) + 2.f * (wdte(0,
    1) + wdte(0, 2)) * (wdte(0, 4) + wdte(0, 5)) + 2.f * wdte(0, 4) * wdte(0,
    5)) / fem::pow2(wdtp(0));
  wids(25, 2) = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) / wdtp(0);
  wids(25, 3) = 0.f;
  FEM_DO_SAFE(i, 0, 40) {
    widp(25, i) = hfac * wdtp(i);
    wide(25, i) = hfac * wdte(i, 0);
  }
  //C
  //C...Z'0:
  zpmas = pmas(32, 1);
  zpfac = aem / (48.f * xw * (1.f - xw)) * zpmas;
  pywidt(cmn, 32, zpmas, wdtp, wdte);
  wids(32, 1) = (fem::pow2((wdte(0, 1) + wdte(0, 2) + wdte(0, 3))) +
    2.f * (wdte(0, 1) + wdte(0, 2)) * (wdte(0, 4) + wdte(0, 5)) + 2.f * wdte(0,
    4) * wdte(0, 5)) / fem::pow2(wdtp(0));
  wids(32, 2) = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) / wdtp(0);
  wids(32, 3) = 0.f;
  FEM_DO_SAFE(i, 0, 40) {
    widp(32, i) = zpfac * wdtp(i);
    wide(32, i) = zpfac * wdte(i, 0);
  }
  //C
  //C...R:
  rmas = pmas(40, 1);
  rfac = 0.08f * rmas / ((mstp(1) - 1) * (1.f + 6.f * (1.f + ulalps(cmn,
    fem::pow2(rmas)) / paru(1))));
  pywidt(cmn, 40, rmas, wdtp, wdte);
  wids(40, 1) = ((wdte(0, 1) + wdte(0, 2)) * (wdte(0, 1) + wdte(0,
    3)) + (wdte(0, 1) + wdte(0, 2) + wdte(0, 1) + wdte(0, 3)) * (wdte(0,
    4) + wdte(0, 5)) + 2.f * wdte(0, 4) * wdte(0, 5)) / fem::pow2(
    wdtp(0));
  wids(40, 2) = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) / wdtp(0);
  wids(40, 3) = (wdte(0, 1) + wdte(0, 3) + wdte(0, 4)) / wdtp(0);
  FEM_DO_SAFE(i, 0, 40) {
    widp(40, i) = wfac * wdtp(i);
    wide(40, i) = wfac * wdte(i, 0);
  }
  //C
  //C...Q:
  kflqm = 1;
  FEM_DO_SAFE(i, 1, fem::min(8, mdcy(21, 3))) {
    idc = i + mdcy(21, 2) - 1;
    if (mdme(idc, 1) <= 0) {
      goto statement_170;
    }
    kflqm = i;
    statement_170:;
  }
  mint(46) = kflqm;
  kfpr(81, 1) = kflqm;
  kfpr(81, 2) = kflqm;
  kfpr(82, 1) = kflqm;
  kfpr(82, 2) = kflqm;
  //C
  //C...Set resonance widths and branching ratios in JETSET.
  FEM_DO_SAFE(i, 1, 6) {
    if (i <= 3) {
      kc = i + 22;
    }
    if (i == 4) {
      kc = 32;
    }
    if (i == 5) {
      kc = 37;
    }
    if (i == 6) {
      kc = 40;
    }
    pmas(kc, 2) = widp(kc, 0);
    pmas(kc, 3) = fem::min(0.9f * pmas(kc, 1), 10.f * pmas(kc, 2));
    FEM_DO_SAFE(j, 1, mdcy(kc, 3)) {
      idc = j + mdcy(kc, 2) - 1;
      brat(idc) = wide(kc, j) / wide(kc, 0);
    }
  }
  //C
  //C...Special cases in treatment of gamma*/Z0: redefine process name.
  if (mstp(43) == 1) {
    proc(1) = "f + fb -> gamma*";
  }
  else if (mstp(43) == 2) {
    proc(1) = "f + fb -> Z0";
  }
  else if (mstp(43) == 3) {
    proc(1) = "f + fb -> gamma*/Z0";
  }
  //C
  //C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name.
  if (mstp(44) == 1) {
    proc(141) = "f + fb -> gamma*";
  }
  else if (mstp(44) == 2) {
    proc(141) = "f + fb -> Z0";
  }
  else if (mstp(44) == 3) {
    proc(141) = "f + fb -> Z'0";
  }
  else if (mstp(44) == 4) {
    proc(141) = "f + fb -> gamma*/Z0";
  }
  else if (mstp(44) == 5) {
    proc(141) = "f + fb -> gamma*/Z'0";
  }
  else if (mstp(44) == 6) {
    proc(141) = "f + fb -> Z0/Z'0";
  }
  else if (mstp(44) == 7) {
    proc(141) = "f + fb -> gamma*/Z0/Z'0";
  }
  //C
}

struct pyxtot_save
{
  arr<float, 2> bcb;
  arr<float> bcc;
  arr<float, 2> bcs;

  pyxtot_save() :
    bcb(dimension(2, 5), fem::fill0),
    bcc(dimension(3), fem::fill0),
    bcs(dimension(5, 8), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
pyxtot(
  common& cmn)
{
  FEM_CMN_SVE(pyxtot);
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  // COMMON pypars
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  // COMMON pyint1
  arr_cref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  // COMMON pyint5
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  // SAVE
  arr_ref<float, 2> bcb(sve.bcb, dimension(2, 5));
  arr_ref<float> bcc(sve.bcc, dimension(3));
  arr_ref<float, 2> bcs(sve.bcs, dimension(5, 8));
  //
  int i = fem::int0;
  int j = fem::int0;
  if (is_called_first_time) {
    {
      static const float values[] = {
        41.74f, 0.66f, 0.0000f, 337.f, 0.0f, 0.0f, -39.3f, 0.48f,
          41.66f, 0.60f, 0.0000f, 306.f, 0.0f, 0.0f, -34.6f, 0.51f,
          41.36f, 0.63f, 0.0000f, 299.f, 7.3f, 0.5f, -40.4f, 0.47f,
          41.68f, 0.63f, 0.0083f, 330.f, 0.0f, 0.0f, -39.0f, 0.48f,
          41.13f, 0.59f, 0.0074f, 278.f, 10.5f, 0.5f, -41.2f, 0.46f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 5) {
        FEM_DO_SAFE(j, 1, 8) {
          data, bcs(i, j);
        }
      }
    }
    {
      static const float values[] = {
        10.79f, -0.049f, 0.040f, 21.5f, 1.23f, 9.92f, -0.027f,
          0.013f, 18.9f, 1.07f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 2) {
        FEM_DO_SAFE(j, 1, 5) {
          data, bcb(i, j);
        }
      }
    }
    {
      static const float values[] = {
        2.0164346f, -0.5590311f, 0.0376279f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        bcc;
    }
  }
  //C
  //C...Parametrizes total, double diffractive, single diffractive and
  //C...elastic cross-sections for different energies and beams.
  //C
  //C...The following data lines are coefficients needed in the
  //C...Block, Cahn parametrization of total cross-section and nuclear
  //C...slope parameter; see below.
  //C
  //C...Total cross-section and nuclear slope parameter for pp and p-pbar
  int nfit = fem::min(5, fem::max(1, mstp(31)));
  float sigp = bcs(nfit, 1) + bcs(nfit, 2) * (-0.25f * fem::pow2(paru(
    1)) * (1.f - 0.25f * bcs(nfit, 3) * fem::pow2(paru(1))) + (1.f +
    0.5f * bcs(nfit, 3) * fem::pow2(paru(1))) * fem::pow2((fem::log(
    vint(2) / bcs(nfit, 4)))) + bcs(nfit, 3) * fem::pow4((fem::log(vint(
    2) / bcs(nfit, 4))))) / (fem::pow2((1.f - 0.25f * bcs(nfit, 3) *
    fem::pow2(paru(1)))) + 2.f * bcs(nfit, 3) * (1.f + 0.25f * bcs(nfit,
    3) * fem::pow2(paru(1))) * fem::pow2((fem::log(vint(2) / bcs(nfit,
    4)))) + fem::pow2(bcs(nfit, 3)) * fem::pow4((fem::log(vint(2) / bcs(nfit,
    4))))) + bcs(nfit, 5) * fem::pow(vint(2), (bcs(nfit, 6) - 1.f)) *
    fem::sin(0.5f * paru(1) * bcs(nfit, 6));
  float sigm = -bcs(nfit, 7) * fem::pow(vint(2), (bcs(nfit, 8) -
    1.f)) * fem::cos(0.5f * paru(1) * bcs(nfit, 8));
  float refp = bcs(nfit, 2) * paru(1) * fem::log(vint(2) / bcs(nfit,
    4)) / (fem::pow2((1.f - 0.25f * bcs(nfit, 3) * fem::pow2(paru(
    1)))) + 2.f * bcs(nfit, 3) * (1.f + 0.25f * bcs(nfit, 3) *
    fem::pow2(paru(1))) + fem::pow2((fem::log(vint(2) / bcs(nfit,
    4)))) + fem::pow2(bcs(nfit, 3)) * fem::pow4((fem::log(vint(2) / bcs(nfit,
    4))))) - bcs(nfit, 5) * fem::pow(vint(2), (bcs(nfit, 6) - 1.f)) *
    fem::cos(0.5f * paru(1) * bcs(nfit, 6));
  float refm = -bcs(nfit, 7) * fem::pow(vint(2), (bcs(nfit, 8) -
    1.f)) * fem::sin(0.5f * paru(1) * bcs(nfit, 8));
  float sigma = sigp - fem::isign(1, mint(11) * mint(12)) * sigm;
  float rho = (refp - fem::isign(1, mint(11) * mint(12)) * refm) / sigma;
  //C
  //C...Nuclear slope parameter B, curvature C:
  nfit = 1;
  if (mstp(31) >= 4) {
    nfit = 2;
  }
  float bp = bcb(nfit, 1) + bcb(nfit, 2) * fem::log(vint(2)) + bcb(nfit,
    3) * fem::pow2((fem::log(vint(2))));
  float bm = bcb(nfit, 4) + bcb(nfit, 5) * fem::log(vint(2));
  float b = bp - fem::isign(1, mint(11) * mint(12)) * sigm / sigp * (bm - bp);
  vint(121) = b;
  float c = -0.5f * bcc(2) / bcc(3) * (1.f - fem::sqrt(fem::max(0.f,
    1.f + 4.f * bcc(3) / fem::pow2(bcc(2)) * (1.e-03f * vint(1) - bcc(
    1)))));
  vint(122) = c;
  //C
  //C...Elastic scattering cross-section (fixed by sigma-tot, rho and B).
  float sigel = fem::pow2(sigma) * (1.f + fem::pow2(rho)) / (16.f *
    paru(1) * paru(5) * b);
  //C
  //C...Single diffractive scattering cross-section from Goulianos:
  float sigsd = 2.f * 0.68f * (1.f + 36.f / vint(2)) * fem::log(
    0.6f + 0.1f * vint(2));
  //C
  //C...Double diffractive scattering cross-section (essentially fixed by
  //C...sigma-sd and sigma-el).
  float sigdd = fem::pow2(sigsd) / (3.f * sigel);
  //C
  //C...Total non-elastic, non-diffractive cross-section.
  float signd = sigma - sigdd - sigsd - sigel;
  //C
  //C...Rescale for pions.
  if (fem::iabs(mint(11)) == 211 && fem::iabs(mint(12)) == 211) {
    sigma = 4.f / 9.f * sigma;
    sigdd = 4.f / 9.f * sigdd;
    sigsd = 4.f / 9.f * sigsd;
    sigel = 4.f / 9.f * sigel;
    signd = 4.f / 9.f * signd;
  }
  else if (fem::iabs(mint(11)) == 211 || fem::iabs(mint(12)) == 211) {
    sigma = 2.f / 3.f * sigma;
    sigdd = 2.f / 3.f * sigdd;
    sigsd = 2.f / 3.f * sigsd;
    sigel = 2.f / 3.f * sigel;
    signd = 2.f / 3.f * signd;
  }
  //C
  //C...Save cross-sections in common block PYPARA.
  vint(101) = sigma;
  vint(102) = sigel;
  vint(103) = sigsd;
  vint(104) = sigdd;
  vint(106) = signd;
  xsec(95, 1) = signd;
  //C
}

//C
//C***********************************************************************
//C
void
pyklim(
  common& cmn,
  int const& ilim)
{
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_cref<float> ckin(cmn.ckin, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  //
  int isub = fem::int0;
  float sqm3 = fem::float0;
  float sqm4 = fem::float0;
  float tau = fem::float0;
  float rm3 = fem::float0;
  float rm4 = fem::float0;
  float be34 = fem::float0;
  float pthmin = fem::float0;
  float yst = fem::float0;
  float cth = fem::float0;
  float taup = fem::float0;
  float x1 = fem::float0;
  float x2 = fem::float0;
  float xf = fem::float0;
  float pth = fem::float0;
  float y3 = fem::float0;
  float y4 = fem::float0;
  float ylarge = fem::float0;
  float ysmall = fem::float0;
  float etalar = fem::float0;
  float etasma = fem::float0;
  float sth = fem::float0;
  float expet3 = fem::float0;
  float expet4 = fem::float0;
  float eta3 = fem::float0;
  float eta4 = fem::float0;
  float cts3 = fem::float0;
  float cts4 = fem::float0;
  float ctslar = fem::float0;
  float ctssma = fem::float0;
  float taumn0 = fem::float0;
  float taumx0 = fem::float0;
  float taumn1 = fem::float0;
  float taumx1 = fem::float0;
  float tm3 = fem::float0;
  float tm4 = fem::float0;
  float ydcosh = fem::float0;
  float taumn2 = fem::float0;
  float taumx2 = fem::float0;
  float cth2mn = fem::float0;
  float cth2mx = fem::float0;
  float taumn3 = fem::float0;
  float taumx3 = fem::float0;
  float taumn4 = fem::float0;
  float taumx4 = fem::float0;
  float taumn5 = fem::float0;
  float taumx5 = fem::float0;
  float taurt = fem::float0;
  float ystmn0 = fem::float0;
  float ystmx0 = fem::float0;
  float ystmn1 = fem::float0;
  float ystmx1 = fem::float0;
  float ystmn2 = fem::float0;
  float ystmx2 = fem::float0;
  float ystmn3 = fem::float0;
  float ystmx3 = fem::float0;
  float yepmn4 = fem::float0;
  float ystmn4 = fem::float0;
  float yepmx4 = fem::float0;
  float ystmx4 = fem::float0;
  float yepsmn = fem::float0;
  float yepsmx = fem::float0;
  float ydifmn = fem::float0;
  float ydifmx = fem::float0;
  float ystmn5 = fem::float0;
  float ystmx5 = fem::float0;
  float cthlim = fem::float0;
  float rzmn = fem::float0;
  float rzmx = fem::float0;
  float yex3mx = fem::float0;
  float yex4mx = fem::float0;
  float yex3mn = fem::float0;
  float yex4mn = fem::float0;
  float ystmn6 = fem::float0;
  float ystmx6 = fem::float0;
  float ctnmn0 = fem::float0;
  float ctnmx0 = fem::float0;
  float ctpmn0 = fem::float0;
  float ctpmx0 = fem::float0;
  float ctnmn1 = fem::float0;
  float ctnmx1 = fem::float0;
  float ctpmn1 = fem::float0;
  float ctpmx1 = fem::float0;
  float ctnmn2 = fem::float0;
  float ctpmx2 = fem::float0;
  float ctnmx2 = fem::float0;
  float ctpmn2 = fem::float0;
  float ctnmn3 = fem::float0;
  float ctnmx3 = fem::float0;
  float ctpmn3 = fem::float0;
  float ctpmx3 = fem::float0;
  float tapmn0 = fem::float0;
  float tapmx0 = fem::float0;
  float tapmn1 = fem::float0;
  float tapmx1 = fem::float0;
  float st2eff = fem::float0;
  //C
  //C...Checks generated variables against pre-set kinematical limits;
  //C...also calculates limits on variables used in generation.
  //C
  //C...Common kinematical expressions.
  isub = mint(1);
  if (isub == 96) {
    goto statement_110;
  }
  sqm3 = vint(63);
  sqm4 = vint(64);
  if (ilim != 1) {
    tau = vint(21);
    rm3 = sqm3 / (tau * vint(2));
    rm4 = sqm4 / (tau * vint(2));
    be34 = fem::sqrt(fem::pow2((1.f - rm3 - rm4)) - 4.f * rm3 * rm4);
  }
  pthmin = ckin(3);
  if (fem::min(sqm3, sqm4) < fem::pow2(ckin(6))) {
    pthmin = fem::max(ckin(3), ckin(5));
  }
  if (ilim == 0) {
    //C...Check generated values of tau, y*, cos(theta-hat), and tau' against
    //C...pre-set kinematical limits.
    yst = vint(22);
    cth = vint(23);
    taup = vint(26);
    if (iset(isub) <= 2) {
      x1 = fem::sqrt(tau) * fem::exp(yst);
      x2 = fem::sqrt(tau) * fem::exp(-yst);
    }
    else {
      x1 = fem::sqrt(taup) * fem::exp(yst);
      x2 = fem::sqrt(taup) * fem::exp(-yst);
    }
    xf = x1 - x2;
    if (tau * vint(2) < fem::pow2(ckin(1))) {
      mint(51) = 1;
    }
    if (ckin(2) >= 0.f && tau * vint(2) > fem::pow2(ckin(2))) {
      mint(51) = 1;
    }
    if (x1 < ckin(21) || x1 > ckin(22)) {
      mint(51) = 1;
    }
    if (x2 < ckin(23) || x2 > ckin(24)) {
      mint(51) = 1;
    }
    if (xf < ckin(25) || xf > ckin(26)) {
      mint(51) = 1;
    }
    if (yst < ckin(7) || yst > ckin(8)) {
      mint(51) = 1;
    }
    if (iset(isub) == 2 || iset(isub) == 4) {
      pth = 0.5f * be34 * fem::sqrt(tau * vint(2) * (1.f - fem::pow2(cth)));
      y3 = yst + 0.5f * fem::log((1.f + rm3 - rm4 + be34 * cth) / (
        1.f + rm3 - rm4 - be34 * cth));
      y4 = yst + 0.5f * fem::log((1.f + rm4 - rm3 - be34 * cth) / (
        1.f + rm4 - rm3 + be34 * cth));
      ylarge = fem::max(y3, y4);
      ysmall = fem::min(y3, y4);
      etalar = 10.f;
      etasma = -10.f;
      sth = fem::sqrt(1.f - fem::pow2(cth));
      if (sth < 1.e-6f) {
        goto statement_100;
      }
      expet3 = ((1.f + rm3 - rm4) * fem::sinh(yst) + be34 * fem::cosh(
        yst) * cth + fem::sqrt(fem::pow2(((1.f + rm3 - rm4) *
        fem::cosh(yst) + be34 * fem::sinh(yst) * cth)) - 4.f *
        rm3)) / (be34 * sth);
      expet4 = ((1.f - rm3 + rm4) * fem::sinh(yst) - be34 * fem::cosh(
        yst) * cth + fem::sqrt(fem::pow2(((1.f - rm3 + rm4) *
        fem::cosh(yst) - be34 * fem::sinh(yst) * cth)) - 4.f *
        rm4)) / (be34 * sth);
      eta3 = fem::log(fem::min(1.e10f, fem::max(1.e-10f, expet3)));
      eta4 = fem::log(fem::min(1.e10f, fem::max(1.e-10f, expet4)));
      etalar = fem::max(eta3, eta4);
      etasma = fem::min(eta3, eta4);
      statement_100:
      cts3 = ((1.f + rm3 - rm4) * fem::sinh(yst) + be34 * fem::cosh(
        yst) * cth) / fem::sqrt(fem::pow2(((1.f + rm3 - rm4) *
        fem::cosh(yst) + be34 * fem::sinh(yst) * cth)) - 4.f * rm3);
      cts4 = ((1.f - rm3 + rm4) * fem::sinh(yst) - be34 * fem::cosh(
        yst) * cth) / fem::sqrt(fem::pow2(((1.f - rm3 + rm4) *
        fem::cosh(yst) - be34 * fem::sinh(yst) * cth)) - 4.f * rm4);
      ctslar = fem::max(cts3, cts4);
      ctssma = fem::min(cts3, cts4);
      if (pth < pthmin) {
        mint(51) = 1;
      }
      if (ckin(4) >= 0.f && pth > ckin(4)) {
        mint(51) = 1;
      }
      if (ylarge < ckin(9) || ylarge > ckin(10)) {
        mint(51) = 1;
      }
      if (ysmall < ckin(11) || ysmall > ckin(12)) {
        mint(51) = 1;
      }
      if (etalar < ckin(13) || etalar > ckin(14)) {
        mint(51) = 1;
      }
      if (etasma < ckin(15) || etasma > ckin(16)) {
        mint(51) = 1;
      }
      if (ctslar < ckin(17) || ctslar > ckin(18)) {
        mint(51) = 1;
      }
      if (ctssma < ckin(19) || ctssma > ckin(20)) {
        mint(51) = 1;
      }
      if (cth < ckin(27) || cth > ckin(28)) {
        mint(51) = 1;
      }
    }
    if (iset(isub) == 3 || iset(isub) == 4) {
      if (taup * vint(2) < fem::pow2(ckin(31))) {
        mint(51) = 1;
      }
      if (ckin(32) >= 0.f && taup * vint(2) > fem::pow2(ckin(32))) {
        mint(51) = 1;
      }
    }
    //C
  }
  else if (ilim == 1) {
    //C...Calculate limits on tau
    //C...0) due to definition
    taumn0 = 0.f;
    taumx0 = 1.f;
    //C...1) due to limits on subsystem mass
    taumn1 = fem::pow2(ckin(1)) / vint(2);
    taumx1 = 1.f;
    if (ckin(2) >= 0.f) {
      taumx1 = fem::pow2(ckin(2)) / vint(2);
    }
    //C...2) due to limits on pT-hat (and non-overlapping rapidity intervals)
    tm3 = fem::sqrt(sqm3 + fem::pow2(pthmin));
    tm4 = fem::sqrt(sqm4 + fem::pow2(pthmin));
    ydcosh = 1.f;
    if (ckin(9) > ckin(12)) {
      ydcosh = fem::cosh(ckin(9) - ckin(12));
    }
    taumn2 = (fem::pow2(tm3) + 2.f * tm3 * tm4 * ydcosh + fem::pow2(
      tm4)) / vint(2);
    taumx2 = 1.f;
    //C...3) due to limits on pT-hat and cos(theta-hat)
    cth2mn = fem::min(fem::pow2(ckin(27)), fem::pow2(ckin(28)));
    cth2mx = fem::max(fem::pow2(ckin(27)), fem::pow2(ckin(28)));
    taumn3 = 0.f;
    if (ckin(27) * ckin(28) > 0.f) {
      taumn3 = fem::pow2((fem::sqrt(sqm3 + fem::pow2(pthmin) / (1.f -
        cth2mn)) + fem::sqrt(sqm4 + fem::pow2(pthmin) / (1.f -
        cth2mn)))) / vint(2);
    }
    taumx3 = 1.f;
    if (ckin(4) >= 0.f && cth2mx < 1.f) {
      taumx3 = fem::pow2((fem::sqrt(sqm3 + fem::pow2(ckin(4)) / (1.f -
        cth2mx)) + fem::sqrt(sqm4 + fem::pow2(ckin(4)) / (1.f -
        cth2mx)))) / vint(2);
    }
    //C...4) due to limits on x1 and x2
    taumn4 = ckin(21) * ckin(23);
    taumx4 = ckin(22) * ckin(24);
    //C...5) due to limits on xF
    taumn5 = 0.f;
    taumx5 = fem::max(1.f - ckin(25), 1.f + ckin(26));
    vint(11) = fem::max(taumn0, taumn1, taumn2, taumn3, taumn4, taumn5);
    vint(31) = fem::min(taumx0, taumx1, taumx2, taumx3, taumx4, taumx5);
    if (mint(43) == 1 && (iset(isub) == 1 || iset(isub) == 2)) {
      vint(11) = 0.99999f;
      vint(31) = 1.00001f;
    }
    if (vint(31) <= vint(11)) {
      mint(51) = 1;
    }
    //C
  }
  else if (ilim == 2) {
    //C...Calculate limits on y*
    if (iset(isub) == 3 || iset(isub) == 4) {
      tau = vint(26);
    }
    taurt = fem::sqrt(tau);
    //C...0) due to kinematics
    ystmn0 = fem::log(taurt);
    ystmx0 = -ystmn0;
    //C...1) due to explicit limits
    ystmn1 = ckin(7);
    ystmx1 = ckin(8);
    //C...2) due to limits on x1
    ystmn2 = fem::log(fem::max(tau, ckin(21)) / taurt);
    ystmx2 = fem::log(fem::max(tau, ckin(22)) / taurt);
    //C...3) due to limits on x2
    ystmn3 = -fem::log(fem::max(tau, ckin(24)) / taurt);
    ystmx3 = -fem::log(fem::max(tau, ckin(23)) / taurt);
    //C...4) due to limits on xF
    yepmn4 = 0.5f * fem::abs(ckin(25)) / taurt;
    ystmn4 = fem::sign(fem::log(fem::sqrt(1.f + fem::pow2(yepmn4)) + yepmn4),
      ckin(25));
    yepmx4 = 0.5f * fem::abs(ckin(26)) / taurt;
    ystmx4 = fem::sign(fem::log(fem::sqrt(1.f + fem::pow2(yepmx4)) + yepmx4),
      ckin(26));
    //C...5) due to simultaneous limits on y-large and y-small
    yepsmn = (rm3 - rm4) * fem::sinh(ckin(9) - ckin(11));
    yepsmx = (rm3 - rm4) * fem::sinh(ckin(10) - ckin(12));
    ydifmn = fem::abs(fem::log(fem::sqrt(1.f + fem::pow2(yepsmn)) - yepsmn));
    ydifmx = fem::abs(fem::log(fem::sqrt(1.f + fem::pow2(yepsmx)) - yepsmx));
    ystmn5 = 0.5f * (ckin(9) + ckin(11) - ydifmn);
    ystmx5 = 0.5f * (ckin(10) + ckin(12) + ydifmx);
    //C...6) due to simultaneous limits on cos(theta-hat) and y-large or
    //C...   y-small
    cthlim = fem::sqrt(1.f - 4.f * fem::pow2(pthmin) / (be34 * tau * vint(2)));
    rzmn = be34 * fem::max(ckin(27), -cthlim);
    rzmx = be34 * fem::min(ckin(28), cthlim);
    yex3mx = (1.f + rm3 - rm4 + rzmx) / fem::max(1e-10f, 1.f + rm3 -
      rm4 - rzmx);
    yex4mx = (1.f + rm4 - rm3 - rzmn) / fem::max(1e-10f, 1.f + rm4 -
      rm3 + rzmn);
    yex3mn = fem::max(1e-10f, 1.f + rm3 - rm4 + rzmn) / (1.f + rm3 -
      rm4 - rzmn);
    yex4mn = fem::max(1e-10f, 1.f + rm4 - rm3 - rzmx) / (1.f + rm4 -
      rm3 + rzmx);
    ystmn6 = ckin(9) - 0.5f * fem::log(fem::max(yex3mx, yex4mx));
    ystmx6 = ckin(12) - 0.5f * fem::log(fem::min(yex3mn, yex4mn));
    vint(12) = fem::max(ystmn0, ystmn1, ystmn2, ystmn3, ystmn4, ystmn5, ystmn6);
    vint(32) = fem::min(ystmx0, ystmx1, ystmx2, ystmx3, ystmx4, ystmx5, ystmx6);
    if (mint(43) == 1) {
      vint(12) = -0.00001f;
      vint(32) = 0.00001f;
    }
    else if (mint(43) == 2) {
      vint(12) = 0.99999f * ystmx0;
      vint(32) = 1.00001f * ystmx0;
    }
    else if (mint(43) == 3) {
      vint(12) = -1.00001f * ystmx0;
      vint(32) = -0.99999f * ystmx0;
    }
    if (vint(32) <= vint(12)) {
      mint(51) = 1;
    }
    //C
  }
  else if (ilim == 3) {
    //C...Calculate limits on cos(theta-hat)
    yst = vint(22);
    //C...0) due to definition
    ctnmn0 = -1.f;
    ctnmx0 = 0.f;
    ctpmn0 = 0.f;
    ctpmx0 = 1.f;
    //C...1) due to explicit limits
    ctnmn1 = fem::min(0.f, ckin(27));
    ctnmx1 = fem::min(0.f, ckin(28));
    ctpmn1 = fem::max(0.f, ckin(27));
    ctpmx1 = fem::max(0.f, ckin(28));
    //C...2) due to limits on pT-hat
    ctnmn2 = -fem::sqrt(1.f - 4.f * fem::pow2(pthmin) / (fem::pow2(
      be34) * tau * vint(2)));
    ctpmx2 = -ctnmn2;
    ctnmx2 = 0.f;
    ctpmn2 = 0.f;
    if (ckin(4) >= 0.f) {
      ctnmx2 = -fem::sqrt(fem::max(0.f, 1.f - 4.f * fem::pow2(ckin(
        4)) / (fem::pow2(be34) * tau * vint(2))));
      ctpmn2 = -ctnmx2;
    }
    //C...3) due to limits on y-large and y-small
    ctnmn3 = fem::min(0.f, fem::max((1.f + rm3 - rm4) / be34 *
      fem::tanh(ckin(11) - yst), -(1.f - rm3 + rm4) / be34 * fem::tanh(
      ckin(10) - yst)));
    ctnmx3 = fem::min(0.f, (1.f + rm3 - rm4) / be34 * fem::tanh(ckin(
      12) - yst), -(1.f - rm3 + rm4) / be34 * fem::tanh(ckin(9) -
      yst));
    ctpmn3 = fem::max(0.f, (1.f + rm3 - rm4) / be34 * fem::tanh(ckin(9) - yst),
      -(1.f - rm3 + rm4) / be34 * fem::tanh(ckin(12) - yst));
    ctpmx3 = fem::max(0.f, fem::min((1.f + rm3 - rm4) / be34 *
      fem::tanh(ckin(10) - yst), -(1.f - rm3 + rm4) / be34 * fem::tanh(
      ckin(11) - yst)));
    vint(13) = fem::max(ctnmn0, ctnmn1, ctnmn2, ctnmn3);
    vint(33) = fem::min(ctnmx0, ctnmx1, ctnmx2, ctnmx3);
    vint(14) = fem::max(ctpmn0, ctpmn1, ctpmn2, ctpmn3);
    vint(34) = fem::min(ctpmx0, ctpmx1, ctpmx2, ctpmx3);
    if (vint(33) <= vint(13) && vint(34) <= vint(14)) {
      mint(51) = 1;
    }
    //C
  }
  else if (ilim == 4) {
    //C...Calculate limits on tau'
    //C...0) due to kinematics
    tapmn0 = tau;
    tapmx0 = 1.f;
    //C...1) due to explicit limits
    tapmn1 = fem::pow2(ckin(31)) / vint(2);
    tapmx1 = 1.f;
    if (ckin(32) >= 0.f) {
      tapmx1 = fem::pow2(ckin(32)) / vint(2);
    }
    vint(16) = fem::max(tapmn0, tapmn1);
    vint(36) = fem::min(tapmx0, tapmx1);
    if (mint(43) == 1) {
      vint(16) = 0.99999f;
      vint(36) = 1.00001f;
    }
    if (vint(36) <= vint(16)) {
      mint(51) = 1;
    }
    //C
  }
  return;
  //C
  //C...Special case for low-pT and multiple interactions:
  //C...effective kinematical limits for tau, y*, cos(theta-hat).
  statement_110:
  if (ilim == 0) {
  }
  else if (ilim == 1) {
    if (mstp(82) <= 1) {
      vint(11) = 4.f * fem::pow2(parp(81)) / vint(2);
    }
    if (mstp(82) >= 2) {
      vint(11) = fem::pow2(parp(82)) / vint(2);
    }
    vint(31) = 1.f;
  }
  else if (ilim == 2) {
    vint(12) = 0.5f * fem::log(vint(21));
    vint(32) = -vint(12);
  }
  else if (ilim == 3) {
    if (mstp(82) <= 1) {
      st2eff = 4.f * fem::pow2(parp(81)) / (vint(21) * vint(2));
    }
    if (mstp(82) >= 2) {
      st2eff = 0.01f * fem::pow2(parp(82)) / (vint(21) * vint(2));
    }
    vint(13) = -fem::sqrt(fem::max(0.f, 1.f - st2eff));
    vint(33) = 0.f;
    vint(14) = 0.f;
    vint(34) = -vint(13);
  }
  //C
}

//C
//C*********************************************************************
//C
void
pykmap(
  common& cmn,
  int const& ivar,
  int const& mvar,
  float const& vvar)
{
  // COMMON pyint1
  arr_cref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  // COMMON pyint2
  arr_cref<int> iset(cmn.iset, dimension(200));
  //
  //C
  //C...Maps a uniform distribution into a distribution of a kinematical
  //C...variable according to one of the possibilities allowed. It is
  //C...assumed that kinematical limits have been set by a PYKLIM call.
  //C
  //C...Convert VVAR to tau variable.
  int isub = mint(1);
  float taumin = fem::float0;
  float taumax = fem::float0;
  float taure = fem::float0;
  float gamre = fem::float0;
  float tau = fem::float0;
  float ratgen = fem::float0;
  float aupp = fem::float0;
  float alow = fem::float0;
  float ystmin = fem::float0;
  float ystmax = fem::float0;
  float yst = fem::float0;
  float rm34 = fem::float0;
  float rsqm = fem::float0;
  float ctnmin = fem::float0;
  float ctnmax = fem::float0;
  float ctpmin = fem::float0;
  float ctpmax = fem::float0;
  float aneg = fem::float0;
  float apos = fem::float0;
  float vctn = fem::float0;
  float cth = fem::float0;
  float vctp = fem::float0;
  float rmnmin = fem::float0;
  float rmnmax = fem::float0;
  float rmpmin = fem::float0;
  float rmpmax = fem::float0;
  float taupmn = fem::float0;
  float taupmx = fem::float0;
  float taup = fem::float0;
  if (ivar == 1) {
    taumin = vint(11);
    taumax = vint(31);
    if (mvar == 3 || mvar == 4) {
      taure = vint(73);
      gamre = vint(74);
    }
    else if (mvar == 5 || mvar == 6) {
      taure = vint(75);
      gamre = vint(76);
    }
    if (mint(43) == 1 && (iset(isub) == 1 || iset(isub) == 2)) {
      tau = 1.f;
    }
    else if (mvar == 1) {
      tau = taumin * fem::pow((taumax / taumin), vvar);
    }
    else if (mvar == 2) {
      tau = taumax * taumin / (taumin + (taumax - taumin) * vvar);
    }
    else if (mvar == 3 || mvar == 5) {
      ratgen = (taure + taumax) / (taure + taumin) * taumin / taumax;
      tau = taure * taumin / ((taure + taumin) * fem::pow(ratgen,
        vvar) - taumin);
    }
    else {
      aupp = fem::atan((taumax - taure) / gamre);
      alow = fem::atan((taumin - taure) / gamre);
      tau = taure + gamre * fem::tan(alow + (aupp - alow) * vvar);
    }
    vint(21) = fem::min(taumax, fem::max(taumin, tau));
    //C
    //C...Convert VVAR to y* variable.
  }
  else if (ivar == 2) {
    ystmin = vint(12);
    ystmax = vint(32);
    if (mint(43) == 1) {
      yst = 0.f;
    }
    else if (mint(43) == 2) {
      if (iset(isub) <= 2) {
        yst = -0.5f * fem::log(vint(21));
      }
      if (iset(isub) >= 3) {
        yst = -0.5f * fem::log(vint(26));
      }
    }
    else if (mint(43) == 3) {
      if (iset(isub) <= 2) {
        yst = 0.5f * fem::log(vint(21));
      }
      if (iset(isub) >= 3) {
        yst = 0.5f * fem::log(vint(26));
      }
    }
    else if (mvar == 1) {
      yst = ystmin + (ystmax - ystmin) * fem::sqrt(vvar);
    }
    else if (mvar == 2) {
      yst = ystmax - (ystmax - ystmin) * fem::sqrt(1.f - vvar);
    }
    else {
      aupp = fem::atan(fem::exp(ystmax));
      alow = fem::atan(fem::exp(ystmin));
      yst = fem::log(fem::tan(alow + (aupp - alow) * vvar));
    }
    vint(22) = fem::min(ystmax, fem::max(ystmin, yst));
    //C
    //C...Convert VVAR to cos(theta-hat) variable.
  }
  else if (ivar == 3) {
    rm34 = 2.f * vint(63) * vint(64) / fem::pow2((vint(21) * vint(2)));
    rsqm = 1.f + rm34;
    if (2.f * fem::pow2(vint(71)) / (vint(21) * vint(2)) < 0.0001f) {
      rm34 = fem::max(rm34, 2.f * fem::pow2(vint(71)) / (vint(21) * vint(2)));
    }
    ctnmin = vint(13);
    ctnmax = vint(33);
    ctpmin = vint(14);
    ctpmax = vint(34);
    if (mvar == 1) {
      aneg = ctnmax - ctnmin;
      apos = ctpmax - ctpmin;
      if (aneg > 0.f && vvar * (aneg + apos) <= aneg) {
        vctn = vvar * (aneg + apos) / aneg;
        cth = ctnmin + (ctnmax - ctnmin) * vctn;
      }
      else {
        vctp = (vvar * (aneg + apos) - aneg) / apos;
        cth = ctpmin + (ctpmax - ctpmin) * vctp;
      }
    }
    else if (mvar == 2) {
      rmnmin = fem::max(rm34, rsqm - ctnmin);
      rmnmax = fem::max(rm34, rsqm - ctnmax);
      rmpmin = fem::max(rm34, rsqm - ctpmin);
      rmpmax = fem::max(rm34, rsqm - ctpmax);
      aneg = fem::log(rmnmin / rmnmax);
      apos = fem::log(rmpmin / rmpmax);
      if (aneg > 0.f && vvar * (aneg + apos) <= aneg) {
        vctn = vvar * (aneg + apos) / aneg;
        cth = rsqm - rmnmin * fem::pow((rmnmax / rmnmin), vctn);
      }
      else {
        vctp = (vvar * (aneg + apos) - aneg) / apos;
        cth = rsqm - rmpmin * fem::pow((rmpmax / rmpmin), vctp);
      }
    }
    else if (mvar == 3) {
      rmnmin = fem::max(rm34, rsqm + ctnmin);
      rmnmax = fem::max(rm34, rsqm + ctnmax);
      rmpmin = fem::max(rm34, rsqm + ctpmin);
      rmpmax = fem::max(rm34, rsqm + ctpmax);
      aneg = fem::log(rmnmax / rmnmin);
      apos = fem::log(rmpmax / rmpmin);
      if (aneg > 0.f && vvar * (aneg + apos) <= aneg) {
        vctn = vvar * (aneg + apos) / aneg;
        cth = rmnmin * fem::pow((rmnmax / rmnmin), vctn) - rsqm;
      }
      else {
        vctp = (vvar * (aneg + apos) - aneg) / apos;
        cth = rmpmin * fem::pow((rmpmax / rmpmin), vctp) - rsqm;
      }
    }
    else if (mvar == 4) {
      rmnmin = fem::max(rm34, rsqm - ctnmin);
      rmnmax = fem::max(rm34, rsqm - ctnmax);
      rmpmin = fem::max(rm34, rsqm - ctpmin);
      rmpmax = fem::max(rm34, rsqm - ctpmax);
      aneg = 1.f / rmnmax - 1.f / rmnmin;
      apos = 1.f / rmpmax - 1.f / rmpmin;
      if (aneg > 0.f && vvar * (aneg + apos) <= aneg) {
        vctn = vvar * (aneg + apos) / aneg;
        cth = rsqm - 1.f / (1.f / rmnmin + aneg * vctn);
      }
      else {
        vctp = (vvar * (aneg + apos) - aneg) / apos;
        cth = rsqm - 1.f / (1.f / rmpmin + apos * vctp);
      }
    }
    else if (mvar == 5) {
      rmnmin = fem::max(rm34, rsqm + ctnmin);
      rmnmax = fem::max(rm34, rsqm + ctnmax);
      rmpmin = fem::max(rm34, rsqm + ctpmin);
      rmpmax = fem::max(rm34, rsqm + ctpmax);
      aneg = 1.f / rmnmin - 1.f / rmnmax;
      apos = 1.f / rmpmin - 1.f / rmpmax;
      if (aneg > 0.f && vvar * (aneg + apos) <= aneg) {
        vctn = vvar * (aneg + apos) / aneg;
        cth = 1.f / (1.f / rmnmin - aneg * vctn) - rsqm;
      }
      else {
        vctp = (vvar * (aneg + apos) - aneg) / apos;
        cth = 1.f / (1.f / rmpmin - apos * vctp) - rsqm;
      }
    }
    if (cth < 0.f) {
      cth = fem::min(ctnmax, fem::max(ctnmin, cth));
    }
    if (cth > 0.f) {
      cth = fem::min(ctpmax, fem::max(ctpmin, cth));
    }
    vint(23) = cth;
    //C
    //C...Convert VVAR to tau' variable.
  }
  else if (ivar == 4) {
    tau = vint(11);
    taupmn = vint(16);
    taupmx = vint(36);
    if (mint(43) == 1) {
      taup = 1.f;
    }
    else if (mvar == 1) {
      taup = taupmn * fem::pow((taupmx / taupmn), vvar);
    }
    else {
      aupp = fem::pow4((1.f - tau / taupmx));
      alow = fem::pow4((1.f - tau / taupmn));
      taup = tau / (1.f - fem::pow((alow + (aupp - alow) * vvar), 0.25f));
    }
    vint(26) = fem::min(taupmx, fem::max(taupmn, taup));
  }
  //C
}

struct pygamm_save
{
  arr<float> b;

  pygamm_save() :
    b(dimension(8), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
float
pygamm(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(pygamm);
  // SAVE
  arr_ref<float> b(sve.b, dimension(8));
  //
  if (is_called_first_time) {
    static const float values[] = {
      -0.57719165f, 0.98820589f, -0.89705694f, 0.91820686f,
        -0.75670408f, 0.48219939f, -0.19352782f, 0.03586834f
    };
    fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
      b;
  }
  //C
  //C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
  //C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
  //C...(Dover, 1965) 6.1.36.
  //Clin      DATA B/-0.577191652,0.988205891,-0.897056937,0.918206857,
  //Clin     &-0.756704078,0.482199394,-0.193527818,0.035868343/
  //C
  int nx = fem::fint(x);
  float dx = x - nx;
  //C
  return_value = 1.f;
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, 8) {
    return_value += b(i) * fem::pow(dx, i);
  }
  int ix = fem::int0;
  if (x < 1.f) {
    return_value = return_value / x;
  }
  else {
    FEM_DO_SAFE(ix, 1, nx - 1) {
      return_value = (x - ix) * return_value;
    }
  }
  //C
  return return_value;
}

struct pystfe_save
{
  arr<fem::str<5> > chdflm;
  fem::str<40> header;
  int init;

  pystfe_save() :
    chdflm(dimension(9), fem::fill0),
    header(fem::char0),
    init(fem::int0)
  {}
};

//C
//C*********************************************************************
//C
void
pystfe(
  common& cmn,
  int const& /* kf */,
  float const& x,
  float const& q2,
  arr_ref<float> xpq)
{
  FEM_CMN_SVE(pystfe);
  xpq(dim1(-6, 6));
  // COMMON ludat2
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  // COMMON pypars
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  //
  // SAVE
  int& init = sve.init;
  //
  str_arr_ref<1> chdflm(sve.chdflm, dimension(9));
  if (is_called_first_time) {
    {
      static const char* values[] = {
        "UPVAL", "DOVAL", "GLUON", "QBAR ", "UBAR ", "SBAR ",
          "CBAR ", "BBAR ", "TBAR "
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chdflm;
    }
    sve.header = "Tung evolution package has been invoked";
    init = 0;
  }
  //C
  //C...This is a dummy routine, where the user can introduce an interface
  //C...to his own external structure function parametrization.
  //C...Arguments in:
  //C...KF : 2212 for p, 211 for pi+; isospin conjugation for n and charge
  //C...    conjugation for pbar, nbar or pi- is performed by PYSTFU.
  //C...X : x value.
  //C...Q2 : Q^2 value.
  //C...Arguments out:
  //C...XPQ(-6:6) : x * f(x,Q2), with index according to KF code,
  //C...    except that gluon is placed in 0. Thus XPQ(0) = xg,
  //C...    XPQ(1) = xd, XPQ(-1) = xdbar, XPQ(2) = xu, XPQ(-2) = xubar,
  //C...    XPQ(3) = xs, XPQ(-3) = xsbar, XPQ(4) = xc, XPQ(-4) = xcbar,
  //C...    XPQ(5) = xb, XPQ(-5) = xbbar, XPQ(6) = xt, XPQ(-6) = xtbar.
  //C...
  //C...One such interface, to the Diemos, Ferroni, Longo, Martinelli
  //C...proton structure functions, already comes with the package. What
  //C...the user needs here is external files with the three routines
  //C...FXG160, FXG260 and FXG360 of the authors above, plus the
  //C...interpolation routine FINT, which is part of the CERN library
  //C...KERNLIB package. To avoid problems with unresolved external
  //C...references, the external calls are commented in the current
  //C...version. To enable this option, remove the C* at the beginning
  //C...of the relevant lines.
  //C...
  //C...Alternatively, the routine can be used as an interface to the
  //C...structure function evolution program of Tung. This can be achieved
  //C...by removing C* at the beginning of some of the lines below.
  //C
  //C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli.
  //C...Allowed variable range 10 GeV2 < Q2 < 1E8 GeV2, 5E-5 < x < .95.
  float xdflm = fem::float0;
  float q2dflm = fem::float0;
  int j = fem::int0;
  arr_1d<9, float> xfdflm(fem::fill0);
  float cxs = fem::float0;
  int i1 = fem::int0;
  int ihdrn = fem::int0;
  int nu = fem::int0;
  int i2 = fem::int0;
  int i3 = fem::int0;
  float alam = fem::float0;
  float tpms = fem::float0;
  float qini = fem::float0;
  float qmax = fem::float0;
  float xmin = fem::float0;
  float q = fem::float0;
  int i = fem::int0;
  float fixq = fem::float0;
  float xps = fem::float0;
  if (mstp(51) >= 11 && mstp(51) <= 13 && mstp(52) <= 1) {
    xdflm = fem::max(0.51e-4f, x);
    q2dflm = fem::max(10.f, fem::min(1e8f, q2));
    if (mstp(52) == 0) {
      q2dflm = 10.f;
    }
    FEM_DO_SAFE(j, 1, 9) {
      if (mstp(52) == 1 && j == 9) {
        q2dflm = q2dflm * fem::pow2((40.f / pmas(6, 1)));
        q2dflm = fem::max(10.f, fem::min(1e8f, q2));
      }
      xfdflm(j) = 0.f;
      //C...Remove C* on following three lines to enable the DFLM options.
      //C*      IF(MSTP(51).EQ.11) CALL FXG160(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
      //C*      IF(MSTP(51).EQ.12) CALL FXG260(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
      //C*      IF(MSTP(51).EQ.13) CALL FXG360(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
    }
    if (x < 0.51e-4f && fem::abs(parp(51) - 1.f) > 0.01f) {
      cxs = fem::pow((0.51e-4f / x), (parp(51) - 1.f));
      FEM_DO_SAFE(j, 1, 7) {
        xfdflm(j) = xfdflm(j) * cxs;
      }
    }
    xpq(0) = xfdflm(3);
    xpq(1) = xfdflm(2) + xfdflm(5);
    xpq(2) = xfdflm(1) + xfdflm(5);
    xpq(3) = xfdflm(6);
    xpq(4) = xfdflm(7);
    xpq(5) = xfdflm(8);
    xpq(6) = xfdflm(9);
    xpq(-1) = xfdflm(5);
    xpq(-2) = xfdflm(5);
    xpq(-3) = xfdflm(6);
    xpq(-4) = xfdflm(7);
    xpq(-5) = xfdflm(8);
    xpq(-6) = xfdflm(9);
    //C
    //C...Proton structure function evolution from Wu-Ki Tung: parton
    //C...distribution functions incorporating heavy quark mass effects.
    //C...Allowed variable range: PARP(52) < Q < PARP(53); PARP(54) < x < 1.
  }
  else {
    if (init == 0) {
      i1 = 0;
      if (mstp(52) == 4) {
        i1 = 1;
      }
      ihdrn = 1;
      nu = mstp(53);
      i2 = mstp(51);
      if (mstp(51) >= 11) {
        i2 = mstp(51) - 3;
      }
      i3 = 0;
      if (mstp(52) == 3) {
        i3 = 1;
      }
      //C
      //C...Convert to Lambda in CWZ scheme (approximately linear relation).
      alam = 0.75f * parp(1);
      tpms = pmas(6, 1);
      qini = parp(52);
      qmax = parp(53);
      xmin = parp(54);
      //C
      //C...Initialize evolution (perform calculation or read results from
      //C...file).
      //C...Remove C* on following two lines to enable Tung initialization.
      //C*        CALL PDFSET(I1,IHDRN,ALAM,TPMS,QINI,QMAX,XMIN,NU,HEADER,
      //C*   &    I2,I3,IRET,IRR)
      init = 1;
    }
    //C
    //C...Put into output array.
    q = fem::sqrt(q2);
    FEM_DO_SAFE(i, -6, 6) {
      fixq = 0.f;
      //C...Remove C* on following line to enable structure function call.
      //C*      FIXQ=MAX(0.,PDF(10,1,I,X,Q,IR))
      xpq(i) = x * fixq;
    }
    //C
    //C...Change order of u and d quarks from Tung to PYTHIA convention.
    xps = xpq(1);
    xpq(1) = xpq(2);
    xpq(2) = xps;
    xps = xpq(-1);
    xpq(-1) = xpq(-2);
    xpq(-2) = xps;
  }
  //C
}

struct pystfu_save
{
  arr<float, 4> cdo;
  arr<float, 5> cehlq;
  arr<float, 4> cow;
  arr<int, 2> nehlq;

  pystfu_save() :
    cdo(dimension(3, 6, 5, 2), fem::fill0),
    cehlq(dimension(6, 6, 2, 8, 2), fem::fill0),
    cow(dimension(3, 5, 4, 2), fem::fill0),
    nehlq(dimension(8, 2), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
pystfu(
  common& cmn,
  int const& kf,
  float const& x,
  float const& q2,
  arr_ref<float> xpq,
  int const& jbt)
{
  FEM_CMN_SVE(pystfu);
  xpq(dim1(-6, 6));
  common_write write(cmn);
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  //
  arr_ref<float, 4> cdo(sve.cdo, dimension(3, 6, 5, 2));
  arr_ref<float, 5> cehlq(sve.cehlq, dimension(6, 6, 2, 8, 2));
  arr_ref<float, 4> cow(sve.cow, dimension(3, 5, 4, 2));
  arr_ref<int, 2> nehlq(sve.nehlq, dimension(8, 2));
  int ix = fem::int0;
  int it = fem::int0;
  int nx = fem::int0;
  int ip = fem::int0;
  int is = fem::int0;
  if (is_called_first_time) {
    {
      static const int values[] = {
        3, 4, 7, 5, 7, 7, 7, 7, 3, 4, 7, 6, 7, 7, 7, 7
      };
      fem::data_of_type<int>(FEM_VALUES_AND_SIZE),
        nehlq;
    }
    {
      static const float values[] = {
        7.677e-01f, -2.087e-01f, -3.303e-01f, -2.517e-02f,
          -1.570e-02f, -1.000e-04f, -5.326e-01f, -2.661e-01f,
          3.201e-01f, 1.192e-01f, 2.434e-02f, 7.620e-03f, 2.162e-01f,
          1.881e-01f, -8.375e-02f, -6.515e-02f, -1.743e-02f,
          -5.040e-03f, -9.211e-02f, -9.952e-02f, 1.373e-02f,
          2.506e-02f, 8.770e-03f, 2.550e-03f, 3.670e-02f, 4.409e-02f,
          9.600e-04f, -7.960e-03f, -3.420e-03f, -1.050e-03f,
          -1.549e-02f, -2.026e-02f, -3.060e-03f, 2.220e-03f,
          1.240e-03f, 4.100e-04f, 2.395e-01f, 2.905e-01f, 9.778e-02f,
          2.149e-02f, 3.440e-03f, 5.000e-04f, 1.751e-02f,
          -6.090e-03f, -2.687e-02f, -1.916e-02f, -7.970e-03f,
          -2.750e-03f, -5.760e-03f, -5.040e-03f, 1.080e-03f,
          2.490e-03f, 1.530e-03f, 7.500e-04f, 1.740e-03f, 1.960e-03f,
          3.000e-04f, -3.400e-04f, -2.900e-04f, -1.800e-04f,
          -5.300e-04f, -6.400e-04f, -1.700e-04f, 4.000e-05f,
          6.000e-05f, 4.000e-05f, 1.700e-04f, 2.200e-04f, 8.000e-05f,
          1.000e-05f, -1.000e-05f, -1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 1, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        7.237e-01f, -2.189e-01f, -2.995e-01f, -1.909e-02f,
          -1.477e-02f, 2.500e-04f, -5.314e-01f, -2.425e-01f,
          3.283e-01f, 1.119e-01f, 2.223e-02f, 7.070e-03f, 2.289e-01f,
          1.890e-01f, -9.859e-02f, -6.900e-02f, -1.747e-02f,
          -5.080e-03f, -1.041e-01f, -1.084e-01f, 2.108e-02f,
          2.975e-02f, 9.830e-03f, 2.830e-03f, 4.394e-02f, 5.116e-02f,
          -1.410e-03f, -1.055e-02f, -4.230e-03f, -1.270e-03f,
          -1.991e-02f, -2.539e-02f, -2.780e-03f, 3.430e-03f,
          1.720e-03f, 5.500e-04f, 2.410e-01f, 2.884e-01f, 9.369e-02f,
          1.900e-02f, 2.530e-03f, 2.400e-04f, 1.765e-02f,
          -9.220e-03f, -3.037e-02f, -2.085e-02f, -8.440e-03f,
          -2.810e-03f, -6.450e-03f, -5.260e-03f, 1.720e-03f,
          3.110e-03f, 1.830e-03f, 8.700e-04f, 2.120e-03f, 2.320e-03f,
          2.600e-04f, -4.900e-04f, -3.900e-04f, -2.300e-04f,
          -6.900e-04f, -8.200e-04f, -2.000e-04f, 7.000e-05f,
          9.000e-05f, 6.000e-05f, 2.400e-04f, 3.100e-04f, 1.100e-04f,
          0.000e+00f, -2.000e-05f, -2.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 1, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        3.813e-01f, -8.090e-02f, -1.634e-01f, -2.185e-02f,
          -8.430e-03f, -6.200e-04f, -2.948e-01f, -1.435e-01f,
          1.665e-01f, 6.638e-02f, 1.473e-02f, 4.080e-03f, 1.252e-01f,
          1.042e-01f, -4.722e-02f, -3.683e-02f, -1.038e-02f,
          -2.860e-03f, -5.478e-02f, -5.678e-02f, 8.900e-03f,
          1.484e-02f, 5.340e-03f, 1.520e-03f, 2.220e-02f, 2.567e-02f,
          -3.000e-05f, -4.970e-03f, -2.160e-03f, -6.500e-04f,
          -9.530e-03f, -1.204e-02f, -1.510e-03f, 1.510e-03f,
          8.300e-04f, 2.700e-04f, 1.261e-01f, 1.354e-01f, 3.958e-02f,
          8.240e-03f, 1.660e-03f, 4.500e-04f, 3.890e-03f,
          -1.159e-02f, -1.625e-02f, -9.610e-03f, -3.710e-03f,
          -1.260e-03f, -1.910e-03f, -5.600e-04f, 1.590e-03f,
          1.590e-03f, 8.400e-04f, 3.900e-04f, 6.400e-04f, 4.900e-04f,
          -1.500e-04f, -2.900e-04f, -1.800e-04f, -1.000e-04f,
          -2.000e-04f, -1.900e-04f, 0.000e+00f, 6.000e-05f,
          4.000e-05f, 3.000e-05f, 7.000e-05f, 8.000e-05f, 2.000e-05f,
          -1.000e-05f, -1.000e-05f, -1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 2, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        3.578e-01f, -8.622e-02f, -1.480e-01f, -1.840e-02f,
          -7.820e-03f, -4.500e-04f, -2.925e-01f, -1.304e-01f,
          1.696e-01f, 6.243e-02f, 1.353e-02f, 3.750e-03f, 1.318e-01f,
          1.041e-01f, -5.486e-02f, -3.872e-02f, -1.038e-02f,
          -2.850e-03f, -6.162e-02f, -6.143e-02f, 1.303e-02f,
          1.740e-02f, 5.940e-03f, 1.670e-03f, 2.643e-02f, 2.957e-02f,
          -1.490e-03f, -6.450e-03f, -2.630e-03f, -7.700e-04f,
          -1.218e-02f, -1.497e-02f, -1.260e-03f, 2.240e-03f,
          1.120e-03f, 3.500e-04f, 1.263e-01f, 1.334e-01f, 3.732e-02f,
          7.070e-03f, 1.260e-03f, 3.400e-04f, 3.660e-03f,
          -1.357e-02f, -1.795e-02f, -1.031e-02f, -3.880e-03f,
          -1.280e-03f, -2.100e-03f, -3.600e-04f, 2.050e-03f,
          1.920e-03f, 9.800e-04f, 4.400e-04f, 7.700e-04f, 5.400e-04f,
          -2.400e-04f, -3.900e-04f, -2.400e-04f, -1.300e-04f,
          -2.600e-04f, -2.300e-04f, 2.000e-05f, 9.000e-05f,
          6.000e-05f, 4.000e-05f, 9.000e-05f, 1.000e-04f, 2.000e-05f,
          -2.000e-05f, -2.000e-05f, -1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 2, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        6.870e-02f, -6.861e-02f, 2.973e-02f, -5.400e-03f, 3.780e-03f,
          -9.700e-04f, -1.802e-02f, 1.400e-04f, 6.490e-03f,
          -8.540e-03f, 1.220e-03f, -1.750e-03f, -4.650e-03f,
          1.480e-03f, -5.930e-03f, 6.000e-04f, -1.030e-03f,
          -8.000e-05f, 6.440e-03f, 2.570e-03f, 2.830e-03f,
          1.150e-03f, 7.100e-04f, 3.300e-04f, -3.930e-03f,
          -2.540e-03f, -1.160e-03f, -7.700e-04f, -3.600e-04f,
          -1.900e-04f, 2.340e-03f, 1.930e-03f, 5.300e-04f,
          3.700e-04f, 1.600e-04f, 9.000e-05f, 1.014e+00f,
          -1.106e+00f, 3.374e-01f, -7.444e-02f, 8.850e-03f,
          -8.700e-04f, 9.233e-01f, -1.285e+00f, 4.475e-01f,
          -9.786e-02f, 1.419e-02f, -1.120e-03f, 4.888e-02f,
          -1.271e-01f, 8.606e-02f, -2.608e-02f, 4.780e-03f,
          -6.000e-04f, -2.691e-02f, 4.887e-02f, -1.771e-02f,
          1.620e-03f, 2.500e-04f, -6.000e-05f, 7.040e-03f,
          -1.113e-02f, 1.590e-03f, 7.000e-04f, -2.000e-04f,
          0.000e+00f, -1.710e-03f, 2.290e-03f, 3.800e-04f,
          -3.500e-04f, 4.000e-05f, 1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 3, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        1.008e-01f, -7.100e-02f, 1.973e-02f, -5.710e-03f, 2.930e-03f,
          -9.900e-04f, -5.271e-02f, -1.823e-02f, 1.792e-02f,
          -6.580e-03f, 1.750e-03f, -1.550e-03f, 1.220e-02f,
          1.763e-02f, -8.690e-03f, -8.800e-04f, -1.160e-03f,
          -2.100e-04f, -1.190e-03f, -7.180e-03f, 2.360e-03f,
          1.890e-03f, 7.700e-04f, 4.100e-04f, -9.100e-04f,
          2.040e-03f, -3.100e-04f, -1.050e-03f, -4.000e-04f,
          -2.400e-04f, 1.190e-03f, -1.700e-04f, -2.000e-04f,
          4.200e-04f, 1.700e-04f, 1.000e-04f, 1.081e+00f,
          -1.189e+00f, 3.868e-01f, -8.617e-02f, 1.115e-02f,
          -1.180e-03f, 9.917e-01f, -1.396e+00f, 4.998e-01f,
          -1.159e-01f, 1.674e-02f, -1.720e-03f, 5.099e-02f,
          -1.338e-01f, 9.173e-02f, -2.885e-02f, 5.890e-03f,
          -6.500e-04f, -3.178e-02f, 5.703e-02f, -2.070e-02f,
          2.440e-03f, 1.100e-04f, -9.000e-05f, 8.970e-03f,
          -1.392e-02f, 2.050e-03f, 6.500e-04f, -2.300e-04f,
          2.000e-05f, -2.340e-03f, 3.010e-03f, 5.000e-04f,
          -3.900e-04f, 6.000e-05f, 1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 3, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        9.482e-01f, -9.578e-01f, 1.009e-01f, -1.051e-01f, 3.456e-02f,
          -3.054e-02f, -9.627e-01f, 5.379e-01f, 3.368e-01f,
          -9.525e-02f, 1.488e-02f, -2.051e-02f, 4.300e-01f,
          -8.306e-02f, -3.372e-01f, 4.902e-02f, -9.160e-03f,
          1.041e-02f, -1.925e-01f, -1.790e-02f, 2.183e-01f,
          7.490e-03f, 4.140e-03f, -1.860e-03f, 8.183e-02f,
          1.926e-02f, -1.072e-01f, -1.944e-02f, -2.770e-03f,
          -5.200e-04f, -3.884e-02f, -1.234e-02f, 5.410e-02f,
          1.879e-02f, 3.350e-03f, 1.040e-03f, 2.948e+01f,
          -3.902e+01f, 1.464e+01f, -3.335e+00f, 5.054e-01f,
          -5.915e-02f, 2.559e+01f, -3.955e+01f, 1.661e+01f,
          -4.299e+00f, 6.904e-01f, -8.243e-02f, -1.663e+00f,
          1.176e+00f, 1.118e+00f, -7.099e-01f, 1.948e-01f,
          -2.404e-02f, -2.168e-01f, 8.170e-01f, -7.169e-01f,
          1.851e-01f, -1.924e-02f, -3.250e-03f, 2.088e-01f,
          -4.355e-01f, 2.239e-01f, -2.446e-02f, -3.620e-03f,
          1.910e-03f, -9.097e-02f, 1.601e-01f, -5.681e-02f,
          -2.500e-03f, 2.580e-03f, -4.700e-04f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 4, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        2.367e+00f, 4.453e-01f, 3.660e-01f, 9.467e-02f, 1.341e-01f,
          1.661e-02f, -3.170e+00f, -1.795e+00f, 3.313e-02f,
          -2.874e-01f, -9.827e-02f, -7.119e-02f, 1.823e+00f,
          1.457e+00f, -2.465e-01f, 3.739e-02f, 6.090e-03f,
          1.814e-02f, -1.033e+00f, -9.827e-01f, 2.136e-01f,
          1.169e-01f, 5.001e-02f, 1.684e-02f, 5.133e-01f, 5.259e-01f,
          -1.173e-01f, -1.139e-01f, -4.988e-02f, -2.021e-02f,
          -2.881e-01f, -3.145e-01f, 5.667e-02f, 9.161e-02f,
          4.568e-02f, 1.951e-02f, 3.036e+01f, -4.062e+01f,
          1.578e+01f, -3.699e+00f, 6.020e-01f, -7.031e-02f,
          2.700e+01f, -4.167e+01f, 1.770e+01f, -4.804e+00f,
          7.862e-01f, -1.060e-01f, -1.909e+00f, 1.357e+00f,
          1.127e+00f, -7.181e-01f, 2.232e-01f, -2.481e-02f,
          -2.488e-01f, 9.781e-01f, -8.127e-01f, 2.094e-01f,
          -2.997e-02f, -4.710e-03f, 2.506e-01f, -5.427e-01f,
          2.672e-01f, -3.103e-02f, -1.800e-03f, 2.870e-03f,
          -1.128e-01f, 2.087e-01f, -6.972e-02f, -2.480e-03f,
          2.630e-03f, -8.400e-04f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 4, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        4.968e-02f, -4.173e-02f, 2.102e-02f, -3.270e-03f, 3.240e-03f,
          -6.700e-04f, -6.150e-03f, -1.294e-02f, 6.740e-03f,
          -6.890e-03f, 9.000e-04f, -1.510e-03f, -8.580e-03f,
          5.050e-03f, -4.900e-03f, -1.600e-04f, -9.400e-04f,
          -1.500e-04f, 7.840e-03f, 1.510e-03f, 2.220e-03f,
          1.400e-03f, 7.000e-04f, 3.500e-04f, -4.410e-03f,
          -2.220e-03f, -8.900e-04f, -8.500e-04f, -3.600e-04f,
          -2.000e-04f, 2.520e-03f, 1.840e-03f, 4.100e-04f,
          3.900e-04f, 1.600e-04f, 9.000e-05f, 9.235e-01f,
          -1.085e+00f, 3.464e-01f, -7.210e-02f, 9.140e-03f,
          -9.100e-04f, 9.315e-01f, -1.274e+00f, 4.512e-01f,
          -9.775e-02f, 1.380e-02f, -1.310e-03f, 4.739e-02f,
          -1.296e-01f, 8.482e-02f, -2.642e-02f, 4.760e-03f,
          -5.700e-04f, -2.653e-02f, 4.953e-02f, -1.735e-02f,
          1.750e-03f, 2.800e-04f, -6.000e-05f, 6.940e-03f,
          -1.132e-02f, 1.480e-03f, 6.500e-04f, -2.100e-04f,
          0.000e+00f, -1.680e-03f, 2.340e-03f, 4.200e-04f,
          -3.400e-04f, 5.000e-05f, 1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 5, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        6.478e-02f, -4.537e-02f, 1.643e-02f, -3.490e-03f, 2.710e-03f,
          -6.700e-04f, -2.223e-02f, -2.126e-02f, 1.247e-02f,
          -6.290e-03f, 1.120e-03f, -1.440e-03f, -1.340e-03f,
          1.362e-02f, -6.130e-03f, -7.900e-04f, -9.000e-04f,
          -2.000e-04f, 5.080e-03f, -3.610e-03f, 1.700e-03f,
          1.830e-03f, 6.800e-04f, 4.000e-04f, -3.580e-03f,
          6.000e-05f, -2.600e-04f, -1.050e-03f, -3.800e-04f,
          -2.300e-04f, 2.420e-03f, 9.300e-04f, -1.000e-04f,
          4.500e-04f, 1.700e-04f, 1.100e-04f, 9.868e-01f,
          -1.171e+00f, 3.940e-01f, -8.459e-02f, 1.124e-02f,
          -1.250e-03f, 1.001e+00f, -1.383e+00f, 5.044e-01f,
          -1.152e-01f, 1.658e-02f, -1.830e-03f, 4.928e-02f,
          -1.368e-01f, 9.021e-02f, -2.935e-02f, 5.800e-03f,
          -6.600e-04f, -3.133e-02f, 5.785e-02f, -2.023e-02f,
          2.630e-03f, 1.600e-04f, -8.000e-05f, 8.840e-03f,
          -1.416e-02f, 1.900e-03f, 5.800e-04f, -2.500e-04f,
          1.000e-05f, -2.300e-03f, 3.080e-03f, 5.500e-04f,
          -3.700e-04f, 7.000e-05f, 1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 5, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        9.270e-03f, -1.817e-02f, 9.590e-03f, -6.390e-03f, 1.690e-03f,
          -1.540e-03f, 5.710e-03f, -1.188e-02f, 6.090e-03f,
          -4.650e-03f, 1.240e-03f, -1.310e-03f, -3.960e-03f,
          7.100e-03f, -3.590e-03f, 1.840e-03f, -3.900e-04f,
          3.400e-04f, 1.120e-03f, -1.960e-03f, 1.120e-03f,
          -4.800e-04f, 1.000e-04f, -4.000e-05f, 4.000e-05f,
          -3.000e-05f, -1.800e-04f, 9.000e-05f, -5.000e-05f,
          -2.000e-05f, -4.200e-04f, 7.300e-04f, -1.600e-04f,
          5.000e-05f, 5.000e-05f, 5.000e-05f, 8.098e-01f,
          -1.042e+00f, 3.398e-01f, -6.824e-02f, 8.760e-03f,
          -9.000e-04f, 8.961e-01f, -1.217e+00f, 4.339e-01f,
          -9.287e-02f, 1.304e-02f, -1.290e-03f, 3.058e-02f,
          -1.040e-01f, 7.604e-02f, -2.415e-02f, 4.600e-03f,
          -5.000e-04f, -2.451e-02f, 4.432e-02f, -1.651e-02f,
          1.430e-03f, 1.200e-04f, -1.000e-04f, 1.122e-02f,
          -1.457e-02f, 2.680e-03f, 5.800e-04f, -1.200e-04f,
          3.000e-05f, -7.730e-03f, 7.330e-03f, -7.600e-04f,
          -2.400e-04f, 1.000e-05f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 6, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        9.980e-03f, -1.945e-02f, 1.055e-02f, -6.870e-03f, 1.860e-03f,
          -1.560e-03f, 5.700e-03f, -1.203e-02f, 6.250e-03f,
          -4.860e-03f, 1.310e-03f, -1.370e-03f, -4.490e-03f,
          7.990e-03f, -4.170e-03f, 2.050e-03f, -4.400e-04f,
          3.300e-04f, 1.470e-03f, -2.480e-03f, 1.460e-03f,
          -5.700e-04f, 1.200e-04f, -1.000e-05f, -9.000e-05f,
          1.500e-04f, -3.200e-04f, 1.200e-04f, -6.000e-05f,
          -4.000e-05f, -4.200e-04f, 7.600e-04f, -1.400e-04f,
          4.000e-05f, 7.000e-05f, 5.000e-05f, 8.698e-01f,
          -1.131e+00f, 3.836e-01f, -8.111e-02f, 1.048e-02f,
          -1.300e-03f, 9.626e-01f, -1.321e+00f, 4.854e-01f,
          -1.091e-01f, 1.583e-02f, -1.700e-03f, 3.057e-02f,
          -1.088e-01f, 8.022e-02f, -2.676e-02f, 5.590e-03f,
          -5.600e-04f, -2.845e-02f, 5.164e-02f, -1.918e-02f,
          2.210e-03f, -4.000e-05f, -1.500e-04f, 1.311e-02f,
          -1.751e-02f, 3.310e-03f, 5.100e-04f, -1.200e-04f,
          5.000e-05f, -8.590e-03f, 8.380e-03f, -9.200e-04f,
          -2.600e-04f, 1.000e-05f, -1.000e-05f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 6, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        9.010e-03f, -1.401e-02f, 7.150e-03f, -4.130e-03f, 1.260e-03f,
          -1.040e-03f, 6.280e-03f, -9.320e-03f, 4.780e-03f,
          -2.890e-03f, 9.100e-04f, -8.200e-04f, -2.930e-03f,
          4.090e-03f, -1.890e-03f, 7.600e-04f, -2.300e-04f,
          1.400e-04f, 3.900e-04f, -1.200e-03f, 4.400e-04f,
          -2.500e-04f, 2.000e-05f, -2.000e-05f, 2.600e-04f,
          1.400e-04f, -8.000e-05f, 1.000e-04f, 1.000e-05f,
          1.000e-05f, -2.600e-04f, 3.200e-04f, 1.000e-05f,
          -1.000e-05f, 1.000e-05f, -1.000e-05f, 8.029e-01f,
          -1.075e+00f, 3.792e-01f, -7.843e-02f, 1.007e-02f,
          -1.090e-03f, 7.903e-01f, -1.099e+00f, 4.153e-01f,
          -9.301e-02f, 1.317e-02f, -1.410e-03f, -1.704e-02f,
          -1.130e-02f, 2.882e-02f, -1.341e-02f, 3.040e-03f,
          -3.600e-04f, -7.200e-04f, 7.230e-03f, -5.160e-03f,
          1.080e-03f, -5.000e-05f, -4.000e-05f, 3.050e-03f,
          -4.610e-03f, 1.660e-03f, -1.300e-04f, -1.000e-05f,
          1.000e-05f, -4.360e-03f, 5.230e-03f, -1.610e-03f,
          2.000e-04f, -2.000e-05f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 7, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        8.980e-03f, -1.459e-02f, 7.510e-03f, -4.410e-03f, 1.310e-03f,
          -1.070e-03f, 5.970e-03f, -9.440e-03f, 4.800e-03f,
          -3.020e-03f, 9.100e-04f, -8.500e-04f, -3.050e-03f,
          4.440e-03f, -2.100e-03f, 8.500e-04f, -2.400e-04f,
          1.400e-04f, 5.300e-04f, -1.300e-03f, 5.600e-04f,
          -2.700e-04f, 3.000e-05f, -2.000e-05f, 2.000e-04f,
          1.400e-04f, -1.100e-04f, 1.000e-04f, 0.000e+00f,
          0.000e+00f, -2.600e-04f, 3.200e-04f, 0.000e+00f,
          -3.000e-05f, 1.000e-05f, -1.000e-05f, 8.672e-01f,
          -1.174e+00f, 4.265e-01f, -9.252e-02f, 1.244e-02f,
          -1.460e-03f, 8.500e-01f, -1.194e+00f, 4.630e-01f,
          -1.083e-01f, 1.614e-02f, -1.830e-03f, -2.241e-02f,
          -5.630e-03f, 2.815e-02f, -1.425e-02f, 3.520e-03f,
          -4.300e-04f, -7.300e-04f, 8.030e-03f, -5.780e-03f,
          1.380e-03f, -1.300e-04f, -4.000e-05f, 3.460e-03f,
          -5.380e-03f, 1.960e-03f, -2.100e-04f, 1.000e-05f,
          1.000e-05f, -4.850e-03f, 5.950e-03f, -1.890e-03f,
          2.600e-04f, -3.000e-05f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 7, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        4.410e-03f, -7.480e-03f, 3.770e-03f, -2.580e-03f, 7.300e-04f,
          -7.100e-04f, 3.840e-03f, -6.050e-03f, 3.030e-03f,
          -2.030e-03f, 5.800e-04f, -5.900e-04f, -8.800e-04f,
          1.660e-03f, -7.500e-04f, 4.700e-04f, -1.000e-04f,
          1.000e-04f, -8.000e-05f, -1.500e-04f, 1.200e-04f,
          -9.000e-05f, 3.000e-05f, 0.000e+00f, 1.300e-04f,
          -2.200e-04f, -2.000e-05f, -2.000e-05f, -2.000e-05f,
          -2.000e-05f, -7.000e-05f, 1.900e-04f, -4.000e-05f,
          2.000e-05f, 0.000e+00f, 0.000e+00f, 6.623e-01f,
          -9.248e-01f, 3.519e-01f, -7.930e-02f, 1.110e-02f,
          -1.180e-03f, 6.380e-01f, -9.062e-01f, 3.582e-01f,
          -8.479e-02f, 1.265e-02f, -1.390e-03f, -2.581e-02f,
          2.125e-02f, 4.190e-03f, -4.980e-03f, 1.490e-03f,
          -2.100e-04f, 7.100e-04f, 5.300e-04f, -1.270e-03f,
          3.900e-04f, -5.000e-05f, -1.000e-05f, 3.850e-03f,
          -5.060e-03f, 1.860e-03f, -3.500e-04f, 4.000e-05f,
          0.000e+00f, -3.530e-03f, 4.460e-03f, -1.500e-03f,
          2.700e-04f, -3.000e-05f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 8, 1);
          }
        }
      }
    }
    {
      static const float values[] = {
        4.260e-03f, -7.530e-03f, 3.830e-03f, -2.680e-03f, 7.600e-04f,
          -7.300e-04f, 3.640e-03f, -6.050e-03f, 3.030e-03f,
          -2.090e-03f, 5.900e-04f, -6.000e-04f, -9.200e-04f,
          1.710e-03f, -8.200e-04f, 5.000e-04f, -1.200e-04f,
          1.000e-04f, -5.000e-05f, -1.600e-04f, 1.300e-04f,
          -9.000e-05f, 3.000e-05f, 0.000e+00f, 1.300e-04f,
          -2.100e-04f, -1.000e-05f, -2.000e-05f, -2.000e-05f,
          -1.000e-05f, -8.000e-05f, 1.800e-04f, -5.000e-05f,
          2.000e-05f, 0.000e+00f, 0.000e+00f, 7.146e-01f,
          -1.007e+00f, 3.932e-01f, -9.246e-02f, 1.366e-02f,
          -1.540e-03f, 6.856e-01f, -9.828e-01f, 3.977e-01f,
          -9.795e-02f, 1.540e-02f, -1.790e-03f, -3.053e-02f,
          2.758e-02f, 2.150e-03f, -4.880e-03f, 1.640e-03f,
          -2.500e-04f, 9.200e-04f, 4.200e-04f, -1.340e-03f,
          4.600e-04f, -8.000e-05f, -1.000e-05f, 4.230e-03f,
          -5.660e-03f, 2.140e-03f, -4.300e-04f, 6.000e-05f,
          0.000e+00f, -3.890e-03f, 5.000e-03f, -1.740e-03f,
          3.300e-04f, -4.000e-05f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(nx, 1, 2) {
        FEM_DO_SAFE(it, 1, 6) {
          FEM_DO_SAFE(ix, 1, 6) {
            data, cehlq(ix, it, nx, 8, 2);
          }
        }
      }
    }
    {
      static const float values[] = {
        4.190e-01f, 3.460e+00f, 4.400e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, 4.000e-03f, 7.240e-01f, -4.860e+00f,
          0.000e+00f, 0.000e+00f, 0.000e+00f, -7.000e-03f,
          -6.600e-02f, 1.330e+00f, 0.000e+00f, 0.000e+00f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 1, 1);
        }
      }
    }
    {
      static const float values[] = {
        3.740e-01f, 3.330e+00f, 6.030e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, 1.400e-02f, 7.530e-01f, -6.220e+00f,
          0.000e+00f, 0.000e+00f, 0.000e+00f, 0.000e+00f,
          -7.600e-02f, 1.560e+00f, 0.000e+00f, 0.000e+00f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 1, 2);
        }
      }
    }
    {
      static const float values[] = {
        7.630e-01f, 4.000e+00f, 0.000e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, -2.370e-01f, 6.270e-01f, -4.210e-01f,
          0.000e+00f, 0.000e+00f, 0.000e+00f, 2.600e-02f,
          -1.900e-02f, 3.300e-02f, 0.000e+00f, 0.000e+00f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 2, 1);
        }
      }
    }
    {
      static const float values[] = {
        7.610e-01f, 3.830e+00f, 0.000e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, -2.320e-01f, 6.270e-01f, -4.180e-01f,
          0.000e+00f, 0.000e+00f, 0.000e+00f, 2.300e-02f,
          -1.900e-02f, 3.600e-02f, 0.000e+00f, 0.000e+00f, 0.000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 2, 2);
        }
      }
    }
    {
      static const float values[] = {
        1.265e+00f, 0.000e+00f, 8.050e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, -1.132e+00f, -3.720e-01f, 1.590e+00f,
          6.310e+00f, -1.050e+01f, 1.470e+01f, 2.930e-01f,
          -2.900e-02f, -1.530e-01f, -2.730e-01f, -3.170e+00f,
          9.800e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 3, 1);
        }
      }
    }
    {
      static const float values[] = {
        1.670e+00f, 0.000e+00f, 9.150e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, -1.920e+00f, -2.730e-01f, 5.300e-01f,
          1.570e+01f, -1.010e+02f, 2.230e+02f, 5.820e-01f,
          -1.640e-01f, -7.630e-01f, -2.830e+00f, 4.470e+01f,
          -1.170e+02f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 3, 2);
        }
      }
    }
    {
      static const float values[] = {
        0.000e+00f, -3.600e-02f, 6.350e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, 1.350e-01f, -2.220e-01f, 3.260e+00f,
          -3.030e+00f, 1.740e+01f, -1.790e+01f, -7.500e-02f,
          -5.800e-02f, -9.090e-01f, 1.500e+00f, -1.130e+01f,
          1.560e+01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 4, 1);
        }
      }
    }
    {
      static const float values[] = {
        0.000e+00f, -1.200e-01f, 3.510e+00f, 0.000e+00f, 0.000e+00f,
          0.000e+00f, 6.700e-02f, -2.330e-01f, 3.660e+00f,
          -4.740e-01f, 9.500e+00f, -1.660e+01f, -3.100e-02f,
          -2.300e-02f, -4.530e-01f, 3.580e-01f, -5.430e+00f,
          1.550e+01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 4, 2);
        }
      }
    }
    {
      static const float values[] = {
        1.560e+00f, 0.000e+00f, 6.000e+00f, 9.000e+00f, 0.000e+00f,
          0.000e+00f, -1.710e+00f, -9.490e-01f, 1.440e+00f,
          -7.190e+00f, -1.650e+01f, 1.530e+01f, 6.380e-01f,
          3.250e-01f, -1.050e+00f, 2.550e-01f, 1.090e+01f, -1.010e+01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 5, 1);
        }
      }
    }
    {
      static const float values[] = {
        8.790e-01f, 0.000e+00f, 4.000e+00f, 9.000e+00f, 0.000e+00f,
          0.000e+00f, -9.710e-01f, -1.160e+00f, 1.230e+00f,
          -5.640e+00f, -7.540e+00f, -5.960e-01f, 4.340e-01f,
          4.760e-01f, -2.540e-01f, -8.170e-01f, 5.500e+00f, 1.260e-01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 6) {
          data, cdo(ip, is, 5, 2);
        }
      }
    }
    {
      static const float values[] = {
        4.0000e-01f, 7.0000e-01f, 0.0000e+00f, 0.0000e+00f,
          0.0000e+00f, -6.2120e-02f, 6.4780e-01f, 0.0000e+00f,
          0.0000e+00f, 0.0000e+00f, -7.1090e-03f, 1.3350e-02f,
          0.0000e+00f, 0.0000e+00f, 0.0000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 1, 1);
        }
      }
    }
    {
      static const float values[] = {
        4.0000e-01f, 6.2800e-01f, 0.0000e+00f, 0.0000e+00f,
          0.0000e+00f, -5.9090e-02f, 6.4360e-01f, 0.0000e+00f,
          0.0000e+00f, 0.0000e+00f, -6.5240e-03f, 1.4510e-02f,
          0.0000e+00f, 0.0000e+00f, 0.0000e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 1, 2);
        }
      }
    }
    {
      static const float values[] = {
        8.8800e-01f, 0.0000e+00f, 3.1100e+00f, 6.0000e+00f,
          0.0000e+00f, -1.8020e+00f, -1.5760e+00f, -1.3170e-01f,
          2.8010e+00f, -1.7280e+01f, 1.8120e+00f, 1.2000e+00f,
          5.0680e-01f, -1.2160e+01f, 2.0490e+01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 2, 1);
        }
      }
    }
    {
      static const float values[] = {
        7.9400e-01f, 0.0000e+00f, 2.8900e+00f, 6.0000e+00f,
          0.0000e+00f, -9.1440e-01f, -1.2370e+00f, 5.9660e-01f,
          -3.6710e+00f, -8.1910e+00f, 5.9660e-01f, 6.5820e-01f,
          -2.5500e-01f, -2.3040e+00f, 7.7580e+00f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 2, 2);
        }
      }
    }
    {
      static const float values[] = {
        9.0000e-01f, 0.0000e+00f, 5.0000e+00f, 0.0000e+00f,
          0.0000e+00f, -2.4280e-01f, -2.1200e-01f, 8.6730e-01f,
          1.2660e+00f, 2.3820e+00f, 1.3860e-01f, 3.6710e-03f,
          4.7470e-02f, -2.2150e+00f, 3.4820e-01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 3, 1);
        }
      }
    }
    {
      static const float values[] = {
        9.0000e-01f, 0.0000e+00f, 5.0000e+00f, 0.0000e+00f,
          0.0000e+00f, -1.4170e-01f, -1.6970e-01f, -2.4740e+00f,
          -2.5340e+00f, 5.6210e-01f, -1.7400e-01f, -9.6230e-02f,
          1.5750e+00f, 1.3780e+00f, -2.7010e-01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 3, 2);
        }
      }
    }
    {
      static const float values[] = {
        0.0000e+00f, -2.2120e-02f, 2.8940e+00f, 0.0000e+00f,
          0.0000e+00f, 7.9280e-02f, -3.7850e-01f, 9.4330e+00f,
          5.2480e+00f, 8.3880e+00f, -6.1340e-02f, -1.0880e-01f,
          -1.0852e+01f, -7.1870e+00f, -1.1610e+01f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 4, 1);
        }
      }
    }
    {
      static const float values[] = {
        0.0000e+00f, -8.8200e-02f, 1.9240e+00f, 0.0000e+00f,
          0.0000e+00f, 6.2290e-02f, -2.8920e-01f, 2.4240e-01f,
          -4.4630e+00f, -8.3670e-01f, -4.0990e-02f, -1.0820e-01f,
          2.0360e+00f, 5.2090e+00f, -4.8400e-02f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(ip, 1, 3) {
        FEM_DO_SAFE(is, 1, 5) {
          data, cow(ip, is, 4, 2);
        }
      }
    }
  }
  float alam = fem::float0;
  int kfl = fem::int0;
  int kfa = fem::int0;
  int kfe = fem::int0;
  int nset = fem::int0;
  float tmin = fem::float0;
  float tmax = fem::float0;
  float t = fem::float0;
  float vt = fem::float0;
  float vx = fem::float0;
  float cxs = fem::float0;
  arr_1d<6, float> tx(fem::fill0);
  arr_1d<6, float> tt(fem::fill0);
  float xqsum = fem::float0;
  arr_1d<6, float> xq(fem::fill0);
  float sd = fem::float0;
  arr_1d<6, float> ts(fem::fill0);
  float eulbt1 = fem::float0;
  float eulbt2 = fem::float0;
  float xps = fem::float0;
  float atnm = fem::float0;
  float bbr2 = fem::float0;
  float abx = fem::float0;
  float apx = fem::float0;
  float aax = fem::float0;
  float rrx = fem::float0;
  static const char* format_1200 =
    "(' Error: bad value of parameter MSTP(51) in PYSTFU,',' MSTP(51) =',i5)";
  //C
  //C                        *******JBT specifies beam or target of the particle
  //C...Gives proton and pi+ parton structure functions according to a few
  //C...different parametrizations. Note that what is coded is x times the
  //C...probability distribution, i.e. xq(x,Q2) etc.
  //C                        ********COMMON BLOCK FROM HIJING
  //C
  //C...The following data lines are coefficients needed in the
  //C...Eichten, Hinchliffe, Lane, Quigg proton structure function
  //C...parametrizations, see below.
  //C...Powers of 1-x in different cases.
  //C...Expansion coefficients for up valence quark distribution.
  //C...Expansion coefficients for down valence quark distribution.
  //C...Expansion coefficients for up and down sea quark distributions.
  //C...Expansion coefficients for gluon distribution.
  //C...Expansion coefficients for strange sea quark distribution.
  //C...Expansion coefficients for charm sea quark distribution.
  //C...Expansion coefficients for bottom sea quark distribution.
  //C...Expansion coefficients for top sea quark distribution.
  //C
  //C...The following data lines are coefficients needed in the
  //C...Duke, Owens proton structure function parametrizations, see below.
  //C...Expansion coefficients for (up+down) valence quark distribution.
  //C...Expansion coefficients for down valence quark distribution.
  //C...Expansion coefficients for (up+down+strange) sea quark distribution.
  //C...Expansion coefficients for charm sea quark distribution.
  //C...Expansion coefficients for gluon distribution.
  //C
  //C...The following data lines are coefficients needed in the
  //C...Owens pion structure function parametrizations, see below.
  //C...Expansion coefficients for up and down valence quark distributions.
  //C...Expansion coefficients for gluon distribution.
  //C...Expansion coefficients for (up+down+strange) quark sea distribution.
  //C...Expansion coefficients for charm quark sea distribution.
  //C
  //C...Euler's beta function, requires ordinary Gamma function
  //Clin-10/25/02 get rid of argument usage mismatch in PYGAMM():
  //C      EULBT(X,Y)=PYGAMM(X)*PYGAMM(Y)/PYGAMM(X+Y)
  //C
  //C...Reset structure functions, check x and hadron flavour.
  alam = 0.f;
  FEM_DO_SAFE(kfl, -6, 6) {
    xpq(kfl) = 0.f;
  }
  if (x < 0.f || x > 1.f) {
    write(mstu(11),
      "(' Error: x value outside physical range, x =',1p,e12.3)"), x;
    return;
  }
  kfa = fem::iabs(kf);
  if (kfa != 211 && kfa != 2212 && kfa != 2112) {
    write(mstu(11),
      "(' Error: illegal particle code for structure function,',' KF =',i5)"),
      kf;
    return;
  }
  //C
  //C...Call user-supplied structure function. Select proton/neutron/pion.
  if (mstp(51) == 0 || mstp(52) >= 2) {
    kfe = kfa;
    if (kfa == 2112) {
      kfe = 2212;
    }
    pystfe(cmn, kfe, x, q2, xpq);
    goto statement_230;
  }
  if (kfa == 211) {
    goto statement_200;
  }
  //C
  if (mstp(51) == 1 || mstp(51) == 2) {
    //C...Proton structure functions from Eichten, Hinchliffe, Lane, Quigg.
    //C...Allowed variable range: 5 GeV2 < Q2 < 1E8 GeV2; 1E-4 < x < 1
    //C
    //C...Determine set, Lamdba and x and t expansion variables.
    nset = mstp(51);
    if (nset == 1) {
      alam = 0.2f;
    }
    if (nset == 2) {
      alam = 0.29f;
    }
    tmin = fem::log(5.f / fem::pow2(alam));
    tmax = fem::log(1e8f / fem::pow2(alam));
    if (mstp(52) == 0) {
      t = tmin;
    }
    else {
      t = fem::log(q2 / fem::pow2(alam));
    }
    vt = fem::max(-1.f, fem::min(1.f, (2.f * t - tmax - tmin) / (tmax - tmin)));
    nx = 1;
    if (x <= 0.1f) {
      nx = 2;
    }
    if (nx == 1) {
      vx = (2.f * x - 1.1f) / 0.9f;
    }
    if (nx == 2) {
      vx = fem::max(-1.f, (2.f * fem::log(x) + 11.51293f) / 6.90776f);
    }
    cxs = 1.f;
    if (x < 1e-4f && fem::abs(parp(51) - 1.f) > 0.01f) {
      cxs = fem::pow((1e-4f / x), (parp(51) - 1.f));
    }
    //C
    //C...Chebyshev polynomials for x and t expansion.
    tx(1) = 1.f;
    tx(2) = vx;
    tx(3) = 2.f * fem::pow2(vx) - 1.f;
    tx(4) = 4.f * fem::pow3(vx) - 3.f * vx;
    tx(5) = 8.f * fem::pow4(vx) - 8.f * fem::pow2(vx) + 1.f;
    tx(6) = 16.f * fem::pow(vx, 5) - 20.f * fem::pow3(vx) + 5.f * vx;
    tt(1) = 1.f;
    tt(2) = vt;
    tt(3) = 2.f * fem::pow2(vt) - 1.f;
    tt(4) = 4.f * fem::pow3(vt) - 3.f * vt;
    tt(5) = 8.f * fem::pow4(vt) - 8.f * fem::pow2(vt) + 1.f;
    tt(6) = 16.f * fem::pow(vt, 5) - 20.f * fem::pow3(vt) + 5.f * vt;
    //C
    //C...Calculate structure functions.
    FEM_DO_SAFE(kfl, 1, 6) {
      xqsum = 0.f;
      FEM_DO_SAFE(it, 1, 6) {
        FEM_DO_SAFE(ix, 1, 6) {
          xqsum += cehlq(ix, it, nx, kfl, nset) * tx(ix) * tt(it);
        }
      }
      xq(kfl) = xqsum * fem::pow((1.f - x), nehlq(kfl, nset)) * cxs;
    }
    //C
    //C...Put into output array.
    xpq(0) = xq(4);
    xpq(1) = xq(2) + xq(3);
    xpq(2) = xq(1) + xq(3);
    xpq(3) = xq(5);
    xpq(4) = xq(6);
    xpq(-1) = xq(3);
    xpq(-2) = xq(3);
    xpq(-3) = xq(5);
    xpq(-4) = xq(6);
    //C
    //C...Special expansion for bottom (thresh effects).
    if (mstp(54) >= 5) {
      if (nset == 1) {
        tmin = 8.1905f;
      }
      if (nset == 2) {
        tmin = 7.4474f;
      }
      if (t <= tmin) {
        goto statement_140;
      }
      vt = fem::max(-1.f, fem::min(1.f, (2.f * t - tmax - tmin) / (
        tmax - tmin)));
      tt(1) = 1.f;
      tt(2) = vt;
      tt(3) = 2.f * fem::pow2(vt) - 1.f;
      tt(4) = 4.f * fem::pow3(vt) - 3.f * vt;
      tt(5) = 8.f * fem::pow4(vt) - 8.f * fem::pow2(vt) + 1.f;
      tt(6) = 16.f * fem::pow(vt, 5) - 20.f * fem::pow3(vt) + 5.f * vt;
      xqsum = 0.f;
      FEM_DO_SAFE(it, 1, 6) {
        FEM_DO_SAFE(ix, 1, 6) {
          xqsum += cehlq(ix, it, nx, 7, nset) * tx(ix) * tt(it);
        }
      }
      xpq(5) = xqsum * fem::pow((1.f - x), nehlq(7, nset));
      xpq(-5) = xpq(5);
      statement_140:;
    }
    //C
    //C...Special expansion for top (thresh effects).
    if (mstp(54) >= 6) {
      if (nset == 1) {
        tmin = 11.5528f;
      }
      if (nset == 2) {
        tmin = 10.8097f;
      }
      tmin += 2.f * fem::log(pmas(6, 1) / 30.f);
      tmax += 2.f * fem::log(pmas(6, 1) / 30.f);
      if (t <= tmin) {
        goto statement_160;
      }
      vt = fem::max(-1.f, fem::min(1.f, (2.f * t - tmax - tmin) / (
        tmax - tmin)));
      tt(1) = 1.f;
      tt(2) = vt;
      tt(3) = 2.f * fem::pow2(vt) - 1.f;
      tt(4) = 4.f * fem::pow3(vt) - 3.f * vt;
      tt(5) = 8.f * fem::pow4(vt) - 8.f * fem::pow2(vt) + 1.f;
      tt(6) = 16.f * fem::pow(vt, 5) - 20.f * fem::pow3(vt) + 5.f * vt;
      xqsum = 0.f;
      FEM_DO_SAFE(it, 1, 6) {
        FEM_DO_SAFE(ix, 1, 6) {
          xqsum += cehlq(ix, it, nx, 8, nset) * tx(ix) * tt(it);
        }
      }
      xpq(6) = xqsum * fem::pow((1.f - x), nehlq(8, nset));
      xpq(-6) = xpq(6);
      statement_160:;
    }
    //C
  }
  else if (mstp(51) == 3 || mstp(51) == 4) {
    //C...Proton structure functions from Duke, Owens.
    //C...Allowed variable range: 4 GeV2 < Q2 < approx 1E6 GeV2.
    //C
    //C...Determine set, Lambda and s expansion parameter.
    nset = mstp(51) - 2;
    if (nset == 1) {
      alam = 0.2f;
    }
    if (nset == 2) {
      alam = 0.4f;
    }
    if (mstp(52) <= 0) {
      sd = 0.f;
    }
    else {
      sd = fem::log(fem::log(fem::max(q2, 4.f) / fem::pow2(alam)) /
        fem::log(4.f / fem::pow2(alam)));
    }
    //C
    //C...Calculate structure functions.
    FEM_DO_SAFE(kfl, 1, 5) {
      FEM_DO_SAFE(is, 1, 6) {
        ts(is) = cdo(1, is, kfl, nset) + cdo(2, is, kfl, nset) * sd + cdo(3,
          is, kfl, nset) * fem::pow2(sd);
      }
      if (kfl <= 2) {
        //C
        //Clin-10/25/02 evaluate EULBT(TS(1),TS(2)+1.):
        //C          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)*(1.+TS(3)*X)/(EULBT(TS(1),
        //C     &    TS(2)+1.)*(1.+TS(3)*TS(1)/(TS(1)+TS(2)+1.)))
        eulbt1 = pygamm(cmn, ts(1)) * pygamm(cmn, ts(2) + 1.f) / pygamm(cmn,
          ts(1) + ts(2) + 1.f);
        xq(kfl) = fem::pow(x, ts(1)) * fem::pow((1.f - x), ts(2)) * (
          1.f + ts(3) * x) / (eulbt1 * (1.f + ts(3) * ts(1) / (ts(
          1) + ts(2) + 1.f)));
      }
      else {
        xq(kfl) = ts(1) * fem::pow(x, ts(2)) * fem::pow((1.f - x), ts(
          3)) * (1.f + ts(4) * x + ts(5) * fem::pow2(x) + ts(6) *
          fem::pow3(x));
      }
      //C
    }
    //C
    //C...Put into output arrays.
    xpq(0) = xq(5);
    xpq(1) = xq(2) + xq(3) / 6.f;
    xpq(2) = 3.f * xq(1) - xq(2) + xq(3) / 6.f;
    xpq(3) = xq(3) / 6.f;
    xpq(4) = xq(4);
    xpq(-1) = xq(3) / 6.f;
    xpq(-2) = xq(3) / 6.f;
    xpq(-3) = xq(3) / 6.f;
    xpq(-4) = xq(4);
    //C
    //C...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli.
    //C...These are accessed via PYSTFE since the files needed may not always
    //C...available.
  }
  else if (mstp(51) >= 11 && mstp(51) <= 13) {
    pystfe(cmn, 2212, x, q2, xpq);
    //C
    //C...Unknown proton parametrization.
  }
  else {
    write(mstu(11), format_1200), mstp(51);
  }
  goto statement_230;
  //C
  statement_200:
  if ((mstp(51) >= 1 && mstp(51) <= 4) || (mstp(51) >= 11 && mstp(51) <= 13)) {
    //C...Pion structure functions from Owens.
    //C...Allowed variable range: 4 GeV2 < Q2 < approx 2000 GeV2.
    //C
    //C...Determine set, Lambda and s expansion variable.
    nset = 1;
    if (mstp(51) == 2 || mstp(51) == 4 || mstp(51) == 13) {
      nset = 2;
    }
    if (nset == 1) {
      alam = 0.2f;
    }
    if (nset == 2) {
      alam = 0.4f;
    }
    if (mstp(52) <= 0) {
      sd = 0.f;
    }
    else {
      sd = fem::log(fem::log(fem::max(q2, 4.f) / fem::pow2(alam)) /
        fem::log(4.f / fem::pow2(alam)));
    }
    //C
    //C...Calculate structure functions.
    FEM_DO_SAFE(kfl, 1, 4) {
      FEM_DO_SAFE(is, 1, 5) {
        ts(is) = cow(1, is, kfl, nset) + cow(2, is, kfl, nset) * sd + cow(3,
          is, kfl, nset) * fem::pow2(sd);
      }
      if (kfl == 1) {
        //C
        //Clin-10/25/02 get rid of argument usage mismatch in PYGAMM():
        //C          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/EULBT(TS(1),TS(2)+1.)
        eulbt2 = pygamm(cmn, ts(1)) * pygamm(cmn, ts(2) + 1.f) / pygamm(cmn,
          ts(1) + ts(2) + 1.f);
        xq(kfl) = fem::pow(x, ts(1)) * fem::pow((1.f - x), ts(2)) / eulbt2;
      }
      else {
        xq(kfl) = ts(1) * fem::pow(x, ts(2)) * fem::pow((1.f - x), ts(
          3)) * (1.f + ts(4) * x + ts(5) * fem::pow2(x));
      }
    }
    //C
    //C...Put into output arrays.
    xpq(0) = xq(2);
    xpq(1) = xq(3) / 6.f;
    xpq(2) = xq(1) + xq(3) / 6.f;
    xpq(3) = xq(3) / 6.f;
    xpq(4) = xq(4);
    xpq(-1) = xq(1) + xq(3) / 6.f;
    xpq(-2) = xq(3) / 6.f;
    xpq(-3) = xq(3) / 6.f;
    xpq(-4) = xq(4);
    //C
    //C...Unknown pion parametrization.
  }
  else {
    write(mstu(11), format_1200), mstp(51);
  }
  //C
  //C...Isospin conjugation for neutron, charge conjugation for antipart.
  statement_230:
  if (kfa == 2112) {
    xps = xpq(1);
    xpq(1) = xpq(2);
    xpq(2) = xps;
    xps = xpq(-1);
    xpq(-1) = xpq(-2);
    xpq(-2) = xps;
  }
  if (kf < 0) {
    FEM_DO_SAFE(kfl, 1, 4) {
      xps = xpq(kfl);
      xpq(kfl) = xpq(-kfl);
      xpq(-kfl) = xps;
    }
  }
  //C
  //C...Check positivity and reset above maximum allowed flavour.
  FEM_DO_SAFE(kfl, -6, 6) {
    xpq(kfl) = fem::max(0.f, xpq(kfl));
    if (fem::iabs(kfl) > mstp(54)) {
      xpq(kfl) = 0.f;
    }
  }
  //C
  //C...consider nuclear effect on the structure function
  if ((jbt != 1 && jbt != 2) || ihpr2(6) == 0 || ihnt2(16) == 1) {
    goto statement_400;
  }
  atnm = ihnt2(2 * jbt - 1);
  if (atnm <= 1.0f) {
    goto statement_400;
  }
  if (jbt == 1) {
    bbr2 = (fem::pow2(yp(1, ihnt2(11))) + fem::pow2(yp(2, ihnt2(
      11)))) / 1.44f / fem::pow(atnm, 0.66666f);
  }
  else if (jbt == 2) {
    bbr2 = (fem::pow2(yt(1, ihnt2(12))) + fem::pow2(yt(2, ihnt2(
      12)))) / 1.44f / fem::pow(atnm, 0.66666f);
  }
  bbr2 = fem::min(1.0f, bbr2);
  abx = (fem::pow(atnm, 0.33333333f) - 1.0f);
  apx = hipr1(6) * 4.0f / 3.0f * abx * fem::sqrt(1.0f - bbr2);
  aax = 1.192f * fem::pow(fem::alog(atnm), 0.1666666f);
  rrx = aax * (fem::pow3(x) - 1.2f * fem::pow2(x) + 0.21f * x) +
    1.0f - (apx - 1.079f * abx * fem::sqrt(x) / fem::alog(atnm +
    1.0f)) * fem::exp(-fem::pow(x, 2.0f) / 0.01f);
  FEM_DO_SAFE(kfl, -6, 6) {
    xpq(kfl) = xpq(kfl) * rrx;
  }
  //C                        ********consider the nuclear effect on the structure
  //C                                function which also depends on the impact
  //C                                parameter of the nuclear reaction
  //C
  statement_400:;
  //C...Formats for error printouts.
  //C
}

//C
//C***********************************************************************
//C
float
pyw1au(
  common& cmn,
  float const& eps,
  int const& ireim)
{
  float return_value = fem::float0;
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  //
  //C
  //C...Calculates real and imaginary parts of the auxiliary function W1;
  //C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
  //C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987
  //C
  //Clin-8/2014:
  //C      ASINH(X)=LOG(X+SQRT(X**2+1.))
  float x = fem::float0;
  acosh(x) = fem::log(x + fem::sqrt(fem::pow2(x) - 1.f));
  //C
  float w1re = fem::float0;
  float w1im = fem::float0;
  if (eps < 0.f) {
    w1re = 2.f * fem::sqrt(1.f - eps) * asinh(fem::sqrt(-1.f / eps));
    w1im = 0.f;
  }
  else if (eps < 1.f) {
    w1re = 2.f * fem::sqrt(1.f - eps) * acosh(fem::sqrt(1.f / eps));
    w1im = -paru(1) * fem::sqrt(1.f - eps);
  }
  else {
    w1re = 2.f * fem::sqrt(eps - 1.f) * fem::asin(fem::sqrt(1.f / eps));
    w1im = 0.f;
  }
  //C
  if (ireim == 1) {
    return_value = w1re;
  }
  if (ireim == 2) {
    return_value = w1im;
  }
  //C
  return return_value;
}

//C
//C***********************************************************************
//C
float
pyw2au(
  common& cmn,
  float const& eps,
  int const& ireim)
{
  float return_value = fem::float0;
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  //
  //C
  //C...Calculates real and imaginary parts of the auxiliary function W2;
  //C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
  //C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987
  //C
  //Clin-8/2014:
  //C      ASINH(X)=LOG(X+SQRT(X**2+1.))
  float x = fem::float0;
  acosh(x) = fem::log(x + fem::sqrt(fem::pow2(x) - 1.f));
  //C
  float w2re = fem::float0;
  float w2im = fem::float0;
  if (eps < 0.f) {
    w2re = 4.f * fem::pow2((asinh(fem::sqrt(-1.f / eps))));
    w2im = 0.f;
  }
  else if (eps < 1.f) {
    w2re = 4.f * fem::pow2((acosh(fem::sqrt(1.f / eps)))) - fem::pow2(paru(1));
    w2im = -4.f * paru(1) * acosh(fem::sqrt(1.f / eps));
  }
  else {
    w2re = -4.f * fem::pow2((fem::asin(fem::sqrt(1.f / eps))));
    w2im = 0.f;
  }
  //C
  if (ireim == 1) {
    return_value = w2re;
  }
  if (ireim == 2) {
    return_value = w2im;
  }
  //C
  return return_value;
}

struct pyspen_save
{
  arr<float> b;

  pyspen_save() :
    b(dim1(0, 14), fem::fill0)
  {}
};

//C
//C***********************************************************************
//C
float
pyspen(
  common& cmn,
  float const& xrein,
  float const& ximin,
  int const& ireim)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(pyspen);
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  //
  // SAVE
  arr_ref<float> b(sve.b, dim1(0, 14));
  //
  if (is_called_first_time) {
    static const float values[] = {
      1.000000e+00f, -5.000000e-01f, 1.666667e-01f, 0.000000e+00f,
        -3.333333e-02f, 0.000000e+00f, 2.380952e-02f, 0.000000e+00f,
        -3.333333e-02f, 0.000000e+00f, 7.575757e-02f, 0.000000e+00f,
        -2.531135e-01f, 0.000000e+00f, 1.166667e+00f
    };
    fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
      b;
  }
  //C
  //C...Calculates real and imaginary part of Spence function; see
  //C...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.
  //C
  float xre = xrein;
  float xim = ximin;
  if (fem::abs(1.f - xre) < 1.e-6f && fem::abs(xim) < 1.e-6f) {
    if (ireim == 1) {
      return_value = fem::pow2(paru(1)) / 6.f;
    }
    if (ireim == 2) {
      return_value = 0.f;
    }
    return return_value;
  }
  //C
  float xmod = fem::sqrt(fem::pow2(xre) + fem::pow2(xim));
  if (xmod < 1.e-6f) {
    if (ireim == 1) {
      return_value = 0.f;
    }
    if (ireim == 2) {
      return_value = 0.f;
    }
    return return_value;
  }
  //C
  float xarg = fem::sign(fem::acos(xre / xmod), xim);
  float sp0re = 0.f;
  float sp0im = 0.f;
  float sgn = 1.f;
  float algxre = fem::float0;
  float algxim = fem::float0;
  if (xmod > 1.f) {
    algxre = fem::log(xmod);
    algxim = xarg - fem::sign(paru(1), xarg);
    sp0re = -fem::pow2(paru(1)) / 6.f - (fem::pow2(algxre) -
      fem::pow2(algxim)) / 2.f;
    sp0im = -algxre * algxim;
    sgn = -1.f;
    xmod = 1.f / xmod;
    xarg = -xarg;
    xre = xmod * fem::cos(xarg);
    xim = xmod * fem::sin(xarg);
  }
  float algyre = fem::float0;
  float algyim = fem::float0;
  if (xre > 0.5f) {
    algxre = fem::log(xmod);
    algxim = xarg;
    xre = 1.f - xre;
    xim = -xim;
    xmod = fem::sqrt(fem::pow2(xre) + fem::pow2(xim));
    xarg = fem::sign(fem::acos(xre / xmod), xim);
    algyre = fem::log(xmod);
    algyim = xarg;
    sp0re += sgn * (fem::pow2(paru(1)) / 6.f - (algxre * algyre -
      algxim * algyim));
    sp0im = sp0im - sgn * (algxre * algyim + algxim * algyre);
    sgn = -sgn;
  }
  //C
  xre = 1.f - xre;
  xim = -xim;
  xmod = fem::sqrt(fem::pow2(xre) + fem::pow2(xim));
  xarg = fem::sign(fem::acos(xre / xmod), xim);
  float zre = -fem::log(xmod);
  float zim = -xarg;
  //C
  float spre = 0.f;
  float spim = 0.f;
  float savere = 1.f;
  float saveim = 0.f;
  int i = fem::int0;
  float termre = fem::float0;
  float termim = fem::float0;
  FEM_DO_SAFE(i, 0, 14) {
    termre = (savere * zre - saveim * zim) / fem::ffloat(i + 1);
    termim = (savere * zim + saveim * zre) / fem::ffloat(i + 1);
    savere = termre;
    saveim = termim;
    spre += b(i) * termre;
    spim += b(i) * termim;
  }
  //C
  if (ireim == 1) {
    return_value = sp0re + sgn * spre;
  }
  if (ireim == 2) {
    return_value = sp0im + sgn * spim;
  }
  //C
  return return_value;
}

//C
//C***********************************************************************
//C
float
pyi3au(
  common& cmn,
  float const& be,
  float const& eps,
  int const& ireim)
{
  float return_value = fem::float0;
  // COMMON ludat1
  arr_cref<float> paru(cmn.paru, dimension(200));
  //
  //C
  //C...Calculates real and imaginary parts of the auxiliary function I3;
  //C...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
  //C...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987
  //C
  float ga = fem::float0;
  if (eps < 1.f) {
    ga = 0.5f * (1.f + fem::sqrt(1.f - eps));
  }
  //C
  float f3re = fem::float0;
  float f3im = fem::float0;
  float rsq = fem::float0;
  float rcthe = fem::float0;
  float rsthe = fem::float0;
  float rcphi = fem::float0;
  float rsphi = fem::float0;
  float r = fem::float0;
  float the = fem::float0;
  float phi = fem::float0;
  if (eps < 0.f) {
    f3re = pyspen(cmn, (ga - 1.f) / (ga + be - 1.f), 0.f, 1) - pyspen(cmn,
      ga / (ga + be - 1.f), 0.f, 1) + pyspen(cmn, (be - ga) / be, 0.f,
      1) - pyspen(cmn, (be - ga) / (be - 1.f), 0.f, 1) + (fem::pow2(
      fem::log(be)) - fem::pow2(fem::log(be - 1.f))) / 2.f + fem::log(
      ga) * fem::log((ga + be - 1.f) / be) + fem::log(ga - 1.f) *
      fem::log((be - 1.f) / (ga + be - 1.f));
    f3im = 0.f;
  }
  else if (eps < 1.f) {
    f3re = pyspen(cmn, (ga - 1.f) / (ga + be - 1.f), 0.f, 1) - pyspen(cmn,
      ga / (ga + be - 1.f), 0.f, 1) + pyspen(cmn, ga / (ga - be),
      0.f, 1) - pyspen(cmn, (ga - 1.f) / (ga - be), 0.f, 1) +
      fem::log(ga / (1.f - ga)) * fem::log((ga + be - 1.f) / (be -
      ga));
    f3im = -paru(1) * fem::log((ga + be - 1.f) / (be - ga));
  }
  else {
    rsq = eps / (eps - 1.f + fem::pow2((2.f * be - 1.f)));
    rcthe = rsq * (1.f - 2.f * be / eps);
    rsthe = fem::sqrt(rsq - fem::pow2(rcthe));
    rcphi = rsq * (1.f + 2.f * (be - 1.f) / eps);
    rsphi = fem::sqrt(rsq - fem::pow2(rcphi));
    r = fem::sqrt(rsq);
    the = fem::acos(rcthe / r);
    phi = fem::acos(rcphi / r);
    f3re = pyspen(cmn, rcthe, rsthe, 1) + pyspen(cmn, rcthe, -rsthe,
      1) - pyspen(cmn, rcphi, rsphi, 1) - pyspen(cmn, rcphi, -rsphi,
      1) + (phi - the) * (phi + the - paru(1));
    f3im = pyspen(cmn, rcthe, rsthe, 2) + pyspen(cmn, rcthe, -rsthe,
      2) - pyspen(cmn, rcphi, rsphi, 2) - pyspen(cmn, rcphi, -rsphi,
      2);
  }
  //C
  if (ireim == 1) {
    return_value = 2.f / (2.f * be - 1.f) * f3re;
  }
  if (ireim == 2) {
    return_value = 2.f / (2.f * be - 1.f) * f3im;
  }
  //C
  return return_value;
}

//C
//C***********************************************************************
//C
void
pysigh(
  common& cmn,
  int& nchn,
  float& sigs)
{
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<float, 2> vckm(cmn.vckm, dimension(4, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> kfin(cmn.kfin, dim1(2).dim2(-40, 40));
  arr_cref<float> ckin(cmn.ckin, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<float> pari(cmn.pari, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_cref<int, 2> kfpr(cmn.kfpr, dimension(200, 2));
  arr_cref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_ref<float, 2> xsfx(cmn.xsfx, dim1(2).dim2(-40, 40));
  arr_ref<int, 2> isig(cmn.isig, dimension(1000, 3));
  arr_ref<float> sigh(cmn.sigh, dimension(1000));
  arr_cref<float, 2> wids(cmn.wids, dim1(21, 40).dim2(3));
  arr_cref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  int isub = fem::int0;
  float taumin = fem::float0;
  float ystmin = fem::float0;
  float ctnmin = fem::float0;
  float ctpmin = fem::float0;
  float xt2min = fem::float0;
  float taupmn = fem::float0;
  float tau = fem::float0;
  float yst = fem::float0;
  float cth = fem::float0;
  float xt2 = fem::float0;
  float taup = fem::float0;
  float taumax = fem::float0;
  float ystmax = fem::float0;
  float ctnmax = fem::float0;
  float ctpmax = fem::float0;
  float xt2max = fem::float0;
  float taupmx = fem::float0;
  arr_1d<2, float> x(fem::fill0);
  float sh = fem::float0;
  float sqm3 = fem::float0;
  float sqm4 = fem::float0;
  float rm3 = fem::float0;
  float rm4 = fem::float0;
  float be34 = fem::float0;
  float rpts = fem::float0;
  float be34l = fem::float0;
  float rm34 = fem::float0;
  float rsqm = fem::float0;
  float rthm = fem::float0;
  float th = fem::float0;
  float uh = fem::float0;
  float sqpth = fem::float0;
  float sh2 = fem::float0;
  float th2 = fem::float0;
  float uh2 = fem::float0;
  float q2 = fem::float0;
  float q2sf = fem::float0;
  int i = fem::int0;
  float xsf = fem::float0;
  arr_1d<13, float> xpq(dim1(-6, 6), fem::fill0);
  int kfl = fem::int0;
  float as = fem::float0;
  float fack = fem::float0;
  float faca = fem::float0;
  float q2as = fem::float0;
  float radc = fem::float0;
  int j = fem::int0;
  arr_2d<2, 81, int> kfac(dim1(2).dim2(-40, 40), fem::fill0);
  int min1 = fem::int0;
  int max1 = fem::int0;
  int min2 = fem::int0;
  int max2 = fem::int0;
  int mina = fem::int0;
  int maxa = fem::int0;
  float sqmz = fem::float0;
  float gmmz = fem::float0;
  float sqmw = fem::float0;
  float gmmw = fem::float0;
  float sqmh = fem::float0;
  float gmmh = fem::float0;
  float sqmzp = fem::float0;
  float gmmzp = fem::float0;
  float sqmhc = fem::float0;
  float gmmhc = fem::float0;
  float sqmr = fem::float0;
  float gmmr = fem::float0;
  float aem = fem::float0;
  float xw = fem::float0;
  float comfac = fem::float0;
  float atau0 = fem::float0;
  float atau1 = fem::float0;
  float h1 = fem::float0;
  float taur1 = fem::float0;
  float gamr1 = fem::float0;
  float atau2 = fem::float0;
  float atau3 = fem::float0;
  float taur2 = fem::float0;
  float gamr2 = fem::float0;
  float atau4 = fem::float0;
  float atau5 = fem::float0;
  float ayst0 = fem::float0;
  float ayst1 = fem::float0;
  float ayst2 = fem::float0;
  float ayst3 = fem::float0;
  float h2 = fem::float0;
  float acth0 = fem::float0;
  float acth1 = fem::float0;
  float acth2 = fem::float0;
  float acth3 = fem::float0;
  float acth4 = fem::float0;
  float h3 = fem::float0;
  float ataup0 = fem::float0;
  float ataup1 = fem::float0;
  float h4 = fem::float0;
  float fzw = fem::float0;
  arr_1d<41, float> wdtp(dim1(0, 40), fem::fill0);
  arr_2d<41, 6, float> wdte(dim1(0, 40).dim2(0, 5), fem::fill0);
  float facz = fem::float0;
  float ei = fem::float0;
  float ai = fem::float0;
  float vi = fem::float0;
  float facf = fem::float0;
  float facw = fem::float0;
  int ia = fem::int0;
  int ja = fem::int0;
  int kchw = fem::int0;
  float fach = fem::float0;
  float rmq = fem::float0;
  float ej = fem::float0;
  float aj = fem::float0;
  float vj = fem::float0;
  float facqq1 = fem::float0;
  float facqqb = fem::float0;
  float facqq2 = fem::float0;
  float facgg1 = fem::float0;
  float facgg2 = fem::float0;
  float facgg = fem::float0;
  float faczg = fem::float0;
  float facwg = fem::float0;
  float fckm = fem::float0;
  float facgz = fem::float0;
  float facgw = fem::float0;
  float faczz = fem::float0;
  float faczw = fem::float0;
  float thuh = fem::float0;
  float visav = fem::float0;
  float aisav = fem::float0;
  float fachz = fem::float0;
  float facww = fem::float0;
  float dsigww = fem::float0;
  float fachw = fem::float0;
  float facqg1 = fem::float0;
  float facqg2 = fem::float0;
  int isde = fem::int0;
  float fgq = fem::float0;
  float facgq = fem::float0;
  float fzq = fem::float0;
  float faczq = fem::float0;
  float facwq = fem::float0;
  float facgg3 = fem::float0;
  float be2 = fem::float0;
  float shang = fem::float0;
  float ashre = fem::float0;
  float ashim = fem::float0;
  float thang = fem::float0;
  float athre = fem::float0;
  float athim = fem::float0;
  float uhang = fem::float0;
  float auhre = fem::float0;
  float auhim = fem::float0;
  float avi = fem::float0;
  float avj = fem::float0;
  float cth2 = fem::float0;
  float atwre = fem::float0;
  float atwim = fem::float0;
  float auwre = fem::float0;
  float auwim = fem::float0;
  float a4re = fem::float0;
  float a4im = fem::float0;
  float ep1 = fem::float0;
  float ep2 = fem::float0;
  float aswre = fem::float0;
  float aswim = fem::float0;
  float be4 = fem::float0;
  float cth3 = fem::float0;
  float sgzang = fem::float0;
  float asgre = fem::float0;
  float asgim = fem::float0;
  float aszre = fem::float0;
  float aszim = fem::float0;
  float tgzang = fem::float0;
  float atgre = fem::float0;
  float atgim = fem::float0;
  float atzre = fem::float0;
  float atzim = fem::float0;
  float alssg = fem::float0;
  int mst115 = fem::int0;
  float q2bn = fem::float0;
  float xrepu = fem::float0;
  float frepu = fem::float0;
  float xattr = fem::float0;
  float fattr = fem::float0;
  float fatre = fem::float0;
  float etare = fem::float0;
  float etaim = fem::float0;
  float eps = fem::float0;
  float root = fem::float0;
  float rln = fem::float0;
  float phire = fem::float0;
  float phiim = fem::float0;
  float eta2 = fem::float0;
  float a5stur = fem::float0;
  float a5stui = fem::float0;
  float sqmq = fem::float0;
  float epss = fem::float0;
  float epsh = fem::float0;
  float facgh = fem::float0;
  float a5tsur = fem::float0;
  float a5tsui = fem::float0;
  float epst = fem::float0;
  float facqh = fem::float0;
  float a2stur = fem::float0;
  float a2stui = fem::float0;
  float a2ustr = fem::float0;
  float a2usti = fem::float0;
  float a2tusr = fem::float0;
  float a2tusi = fem::float0;
  float a4stur = fem::float0;
  float a4stui = fem::float0;
  float epsu = fem::float0;
  float bestu = fem::float0;
  float beust = fem::float0;
  float betus = fem::float0;
  float beuts = fem::float0;
  float betsu = fem::float0;
  float besut = fem::float0;
  float w3stur = fem::float0;
  float w3stui = fem::float0;
  float w3sutr = fem::float0;
  float w3suti = fem::float0;
  float w3tsur = fem::float0;
  float w3tsui = fem::float0;
  float w3tusr = fem::float0;
  float w3tusi = fem::float0;
  float w3ustr = fem::float0;
  float w3usti = fem::float0;
  float w3utsr = fem::float0;
  float w3utsi = fem::float0;
  float b2stur = fem::float0;
  float b2stui = fem::float0;
  float b2sutr = fem::float0;
  float b2suti = fem::float0;
  float b2tsur = fem::float0;
  float b2tsui = fem::float0;
  float b2tusr = fem::float0;
  float b2tusi = fem::float0;
  float b2ustr = fem::float0;
  float b2usti = fem::float0;
  float b2utsr = fem::float0;
  float b2utsi = fem::float0;
  float b4stur = fem::float0;
  float b4stui = fem::float0;
  float b4tusr = fem::float0;
  float b4tusi = fem::float0;
  float b4ustr = fem::float0;
  float b4usti = fem::float0;
  float asre = fem::float0;
  float asim = fem::float0;
  float a0stur = fem::float0;
  float a0stui = fem::float0;
  float a0tsur = fem::float0;
  float a0tsui = fem::float0;
  float a0utsr = fem::float0;
  float a0utsi = fem::float0;
  float a1stur = fem::float0;
  float a1stui = fem::float0;
  float faczp = fem::float0;
  float api = fem::float0;
  float vpi = fem::float0;
  float fhc = fem::float0;
  int il = fem::int0;
  int iu = fem::int0;
  float rmql = fem::float0;
  float rmqu = fem::float0;
  float fachc = fem::float0;
  int kchhc = fem::int0;
  float facr = fem::float0;
  float fhcq = fem::float0;
  float fachcq = fem::float0;
  int ichn = fem::int0;
  int kfl1 = fem::int0;
  int kfl2 = fem::int0;
  //C
  //C...Differential matrix elements for all included subprocesses.
  //C...Note that what is coded is (disregarding the COMFAC factor)
  //C...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,
  //C...when d(sigma-hat) is given in the zero-width limit, the delta
  //C...function in tau is replaced by a Breit-Wigner:
  //C...1/pi*(s*m_res*Gamma_res)/((s*tau-m_res^2)^2+(m_res*Gamma_res)^2);
  //C...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);
  //C...i.e., dimensionless quantities. COMFAC contains the factor
  //C...pi/s and the conversion factor from GeV^-2 to mb.
  //C
  //C...Reset number of channels and cross-section.
  nchn = 0;
  sigs = 0.f;
  //C
  //C...Read kinematical variables and limits.
  isub = mint(1);
  taumin = vint(11);
  ystmin = vint(12);
  ctnmin = vint(13);
  ctpmin = vint(14);
  xt2min = vint(15);
  taupmn = vint(16);
  tau = vint(21);
  yst = vint(22);
  cth = vint(23);
  xt2 = vint(25);
  taup = vint(26);
  taumax = vint(31);
  ystmax = vint(32);
  ctnmax = vint(33);
  ctpmax = vint(34);
  xt2max = vint(35);
  taupmx = vint(36);
  //C
  //C...Derive kinematical quantities.
  if (iset(isub) <= 2 || iset(isub) == 5) {
    x(1) = fem::sqrt(tau) * fem::exp(yst);
    x(2) = fem::sqrt(tau) * fem::exp(-yst);
  }
  else {
    x(1) = fem::sqrt(taup) * fem::exp(yst);
    x(2) = fem::sqrt(taup) * fem::exp(-yst);
  }
  if (mint(43) == 4 && iset(isub) >= 1 && (x(1) > 0.999f || x(2) > 0.999f)) {
    return;
  }
  sh = tau * vint(2);
  sqm3 = vint(63);
  sqm4 = vint(64);
  rm3 = sqm3 / sh;
  rm4 = sqm4 / sh;
  be34 = fem::sqrt(fem::pow2((1.f - rm3 - rm4)) - 4.f * rm3 * rm4);
  rpts = 4.f * fem::pow2(vint(71)) / sh;
  be34l = fem::sqrt(fem::max(0.f, fem::pow2((1.f - rm3 - rm4)) -
    4.f * rm3 * rm4 - rpts));
  rm34 = 2.f * rm3 * rm4;
  rsqm = 1.f + rm34;
  rthm = (4.f * rm3 * rm4 + rpts) / (1.f - rm3 - rm4 + be34l);
  th = -0.5f * sh * fem::max(rthm, 1.f - rm3 - rm4 - be34 * cth);
  uh = -0.5f * sh * fem::max(rthm, 1.f - rm3 - rm4 + be34 * cth);
  sqpth = 0.25f * sh * fem::pow2(be34) * (1.f - fem::pow2(cth));
  sh2 = fem::pow2(sh);
  th2 = fem::pow2(th);
  uh2 = fem::pow2(uh);
  //C
  //C...Choice of Q2 scale.
  if (iset(isub) == 1 || iset(isub) == 3) {
    q2 = sh;
  }
  else if (fem::mod(iset(isub), 2) == 0 || iset(isub) == 5) {
    if (mstp(32) == 1) {
      q2 = 2.f * sh * th * uh / (fem::pow2(sh) + fem::pow2(th) + fem::pow2(uh));
    }
    else if (mstp(32) == 2) {
      q2 = sqpth + 0.5f * (sqm3 + sqm4);
    }
    else if (mstp(32) == 3) {
      q2 = fem::min(-th, -uh);
    }
    else if (mstp(32) == 4) {
      q2 = sh;
    }
    if (iset(isub) == 5 && mstp(82) >= 2) {
      q2 += fem::pow2(parp(82));
    }
  }
  //C
  //C...Store derived kinematical quantities.
  vint(41) = x(1);
  vint(42) = x(2);
  vint(44) = sh;
  vint(43) = fem::sqrt(sh);
  vint(45) = th;
  vint(46) = uh;
  vint(48) = sqpth;
  vint(47) = fem::sqrt(sqpth);
  vint(50) = taup * vint(2);
  vint(49) = fem::sqrt(fem::max(0.f, vint(50)));
  vint(52) = q2;
  vint(51) = fem::sqrt(q2);
  //C
  //C...Calculate parton structure functions.
  if (iset(isub) <= 0) {
    goto statement_145;
  }
  if (mint(43) >= 2) {
    q2sf = q2;
    if (iset(isub) == 3 || iset(isub) == 4) {
      q2sf = fem::pow2(pmas(23, 1));
      if (isub == 8 || isub == 76 || isub == 77) {
        q2sf = fem::pow2(pmas(24, 1));
      }
    }
    FEM_DO_SAFE(i, 3 - mint(41), mint(42)) {
      xsf = x(i);
      if (iset(isub) == 5) {
        xsf = x(i) / vint(142 + i);
      }
      pystfu(cmn, mint(10 + i), xsf, q2sf, xpq, i);
      FEM_DO_SAFE(kfl, -6, 6) {
        xsfx(i, kfl) = xpq(kfl);
      }
    }
  }
  //C
  //C...Calculate alpha_strong and K-factor.
  if (mstp(33) != 3) {
    as = ulalps(cmn, q2);
  }
  fack = 1.f;
  faca = 1.f;
  if (mstp(33) == 1) {
    fack = parp(31);
  }
  else if (mstp(33) == 2) {
    fack = parp(31);
    faca = parp(32) / parp(31);
  }
  else if (mstp(33) == 3) {
    q2as = parp(33) * q2;
    if (iset(isub) == 5 && mstp(82) >= 2) {
      q2as += paru(112) * parp(82);
    }
    as = ulalps(cmn, q2as);
  }
  radc = 1.f + as / paru(1);
  //C
  //C...Set flags for allowed reacting partons/leptons.
  FEM_DO_SAFE(i, 1, 2) {
    FEM_DO_SAFE(j, -40, 40) {
      kfac(i, j) = 0;
    }
    if (mint(40 + i) == 1) {
      kfac(i, mint(10 + i)) = 1;
    }
    else {
      FEM_DO_SAFE(j, -40, 40) {
        kfac(i, j) = kfin(i, j);
        if (fem::abs(j) > mstp(54) && j != 21) {
          kfac(i, j) = 0;
        }
        if (fem::abs(j) <= 6) {
          if (xsfx(i, j) < 1.e-10f) {
            kfac(i, j) = 0;
          }
        }
        else if (j == 21) {
          if (xsfx(i, 0) < 1.e-10f) {
            kfac(i, 21) = 0;
          }
        }
      }
    }
  }
  //C
  //C...Lower and upper limit for flavour loops.
  min1 = 0;
  max1 = 0;
  min2 = 0;
  max2 = 0;
  FEM_DO_SAFE(j, -20, 20) {
    if (kfac(1, -j) == 1) {
      min1 = -j;
    }
    if (kfac(1, j) == 1) {
      max1 = j;
    }
    if (kfac(2, -j) == 1) {
      min2 = -j;
    }
    if (kfac(2, j) == 1) {
      max2 = j;
    }
  }
  mina = fem::min(min1, min2);
  maxa = fem::max(max1, max2);
  //C
  //C...Common conversion factors (including Jacobian) for subprocesses.
  sqmz = fem::pow2(pmas(23, 1));
  gmmz = pmas(23, 1) * pmas(23, 2);
  sqmw = fem::pow2(pmas(24, 1));
  gmmw = pmas(24, 1) * pmas(24, 2);
  sqmh = fem::pow2(pmas(25, 1));
  gmmh = pmas(25, 1) * pmas(25, 2);
  sqmzp = fem::pow2(pmas(32, 1));
  gmmzp = pmas(32, 1) * pmas(32, 2);
  sqmhc = fem::pow2(pmas(37, 1));
  gmmhc = pmas(37, 1) * pmas(37, 2);
  sqmr = fem::pow2(pmas(40, 1));
  gmmr = pmas(40, 1) * pmas(40, 2);
  aem = paru(101);
  xw = paru(102);
  //C
  //C...Phase space integral in tau and y*.
  comfac = paru(1) * paru(5) / vint(2);
  if (mint(43) == 4) {
    comfac = comfac * fack;
  }
  if ((mint(43) >= 2 || iset(isub) == 3 || iset(isub) == 4) && iset(
      isub) != 5) {
    atau0 = fem::log(taumax / taumin);
    atau1 = (taumax - taumin) / (taumax * taumin);
    h1 = coef(isub, 1) + (atau0 / atau1) * coef(isub, 2) / tau;
    if (mint(72) >= 1) {
      taur1 = vint(73);
      gamr1 = vint(74);
      atau2 = fem::log(taumax / taumin * (taumin + taur1) / (taumax +
        taur1)) / taur1;
      atau3 = (fem::atan((taumax - taur1) / gamr1) - fem::atan((
        taumin - taur1) / gamr1)) / gamr1;
      h1 += (atau0 / atau2) * coef(isub, 3) / (tau + taur1) + (atau0 /
        atau3) * coef(isub, 4) * tau / (fem::pow2((tau - taur1)) +
        fem::pow2(gamr1));
    }
    if (mint(72) == 2) {
      taur2 = vint(75);
      gamr2 = vint(76);
      atau4 = fem::log(taumax / taumin * (taumin + taur2) / (taumax +
        taur2)) / taur2;
      atau5 = (fem::atan((taumax - taur2) / gamr2) - fem::atan((
        taumin - taur2) / gamr2)) / gamr2;
      h1 += (atau0 / atau4) * coef(isub, 5) / (tau + taur2) + (atau0 /
        atau5) * coef(isub, 6) * tau / (fem::pow2((tau - taur2)) +
        fem::pow2(gamr2));
    }
    comfac = comfac * atau0 / (tau * h1);
  }
  if (mint(43) == 4 && iset(isub) != 5) {
    ayst0 = ystmax - ystmin;
    ayst1 = 0.5f * fem::pow2((ystmax - ystmin));
    ayst2 = ayst1;
    ayst3 = 2.f * (fem::atan(fem::exp(ystmax)) - fem::atan(fem::exp(ystmin)));
    h2 = (ayst0 / ayst1) * coef(isub, 7) * (yst - ystmin) + (ayst0 /
      ayst2) * coef(isub, 8) * (ystmax - yst) + (ayst0 / ayst3) * coef(isub,
      9) / fem::cosh(yst);
    comfac = comfac * ayst0 / h2;
  }
  //C
  //C...2 -> 1 processes: reduction in angular part of phase space integral
  //C...for case of decaying resonance.
  acth0 = ctnmax - ctnmin + ctpmax - ctpmin;
  //Clin-4/2008 modified a la pythia6115.f to avoid invalid MDCY subcript#1,
  //C     also break up compound IF statements:
  //C      IF((ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3).AND.
  //C     &MDCY(KFPR(ISUB,1),1).EQ.1) THEN
  //C        IF(KFPR(ISUB,1).EQ.25.OR.KFPR(ISUB,1).EQ.37) THEN
  //C          COMFAC=COMFAC*0.5*ACTH0
  //C        ELSE
  //C          COMFAC=COMFAC*0.125*(3.*ACTH0+CTNMAX**3-CTNMIN**3+
  //C     &    CTPMAX**3-CTPMIN**3)
  //C        ENDIF
  if (iset(isub) == 1 || iset(isub) == 3) {
    if (mdcy(lucomp(cmn, kfpr(isub, 1)), 1) == 1) {
      if (kfpr(isub, 1) == 25 || kfpr(isub, 1) == 37) {
        comfac = comfac * 0.5f * acth0;
      }
      else {
        comfac = comfac * 0.125f * (3.f * acth0 + fem::pow3(ctnmax) -
          fem::pow3(ctnmin) + fem::pow3(ctpmax) - fem::pow3(ctpmin));
      }
    }
    //C
    //C...2 -> 2 processes: angular part of phase space integral.
  }
  else if (iset(isub) == 2 || iset(isub) == 4) {
    acth1 = fem::log((fem::max(rm34, rsqm - ctnmin) * fem::max(rm34,
      rsqm - ctpmin)) / (fem::max(rm34, rsqm - ctnmax) * fem::max(rm34,
      rsqm - ctpmax)));
    acth2 = fem::log((fem::max(rm34, rsqm + ctnmax) * fem::max(rm34,
      rsqm + ctpmax)) / (fem::max(rm34, rsqm + ctnmin) * fem::max(rm34,
      rsqm + ctpmin)));
    acth3 = 1.f / fem::max(rm34, rsqm - ctnmax) - 1.f / fem::max(rm34,
      rsqm - ctnmin) + 1.f / fem::max(rm34, rsqm - ctpmax) - 1.f /
      fem::max(rm34, rsqm - ctpmin);
    acth4 = 1.f / fem::max(rm34, rsqm + ctnmin) - 1.f / fem::max(rm34,
      rsqm + ctnmax) + 1.f / fem::max(rm34, rsqm + ctpmin) - 1.f /
      fem::max(rm34, rsqm + ctpmax);
    h3 = coef(isub, 10) + (acth0 / acth1) * coef(isub, 11) / fem::max(rm34,
      rsqm - cth) + (acth0 / acth2) * coef(isub, 12) / fem::max(rm34,
      rsqm + cth) + (acth0 / acth3) * coef(isub, 13) / fem::pow2(fem::max(rm34,
      rsqm - cth)) + (acth0 / acth4) * coef(isub, 14) / fem::pow2(
      fem::max(rm34, rsqm + cth));
    comfac = comfac * acth0 * 0.5f * be34 / h3;
  }
  //C
  //C...2 -> 3, 4 processes: phace space integral in tau'.
  if (mint(43) >= 2 && (iset(isub) == 3 || iset(isub) == 4)) {
    ataup0 = fem::log(taupmx / taupmn);
    ataup1 = (fem::pow4((1.f - tau / taupmx)) - fem::pow4((1.f -
      tau / taupmn))) / (4.f * tau);
    h4 = coef(isub, 15) + ataup0 / ataup1 * coef(isub, 16) / taup *
      fem::pow3((1.f - tau / taup));
    if (1.f - tau / taup > 1.e-4f) {
      fzw = (1.f + tau / taup) * fem::log(taup / tau) - 2.f * (1.f -
        tau / taup);
    }
    else {
      fzw = 1.f / 6.f * fem::pow3((1.f - tau / taup)) * tau / taup;
    }
    comfac = comfac * ataup0 * fzw / h4;
  }
  //C
  //C...Phase space integral for low-pT and multiple interactions.
  if (iset(isub) == 5) {
    comfac = paru(1) * paru(5) * fack * 0.5f * vint(2) / sh2;
    atau0 = fem::log(2.f * (1.f + fem::sqrt(1.f - xt2)) / xt2 - 1.f);
    atau1 = 2.f * fem::atan(1.f / xt2 - 1.f) / fem::sqrt(xt2);
    h1 = coef(isub, 1) + (atau0 / atau1) * coef(isub, 2) / fem::sqrt(tau);
    comfac = comfac * atau0 / h1;
    ayst0 = ystmax - ystmin;
    ayst1 = 0.5f * fem::pow2((ystmax - ystmin));
    ayst3 = 2.f * (fem::atan(fem::exp(ystmax)) - fem::atan(fem::exp(ystmin)));
    h2 = (ayst0 / ayst1) * coef(isub, 7) * (yst - ystmin) + (ayst0 /
      ayst1) * coef(isub, 8) * (ystmax - yst) + (ayst0 / ayst3) * coef(isub,
      9) / fem::cosh(yst);
    comfac = comfac * ayst0 / h2;
    if (mstp(82) <= 1) {
      comfac = comfac * fem::pow2(xt2) * (1.f / vint(149) - 1.f);
    }
    //C...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is
    //C...introduced to make cross-section finite for xT2 -> 0.
    if (mstp(82) >= 2) {
      comfac = comfac * fem::pow2(xt2) / (vint(149) * (1.f + vint(149)));
    }
  }
  //C
  //C...A: 2 -> 1, tree diagrams.
  //C
  statement_145:
  if (isub <= 10) {
    if (isub == 1) {
      //C...f + fb -> gamma*/Z0.
      mint(61) = 2;
      pywidt(cmn, 23, fem::sqrt(sh), wdtp, wdte);
      facz = comfac * fem::pow2(aem) * 4.f / 3.f;
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_150;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        facf = 1.f;
        if (fem::iabs(i) <= 10) {
          facf = faca / 3.f;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facf * facz * (fem::pow2(ei) * vint(111) + ei *
          vi / (8.f * xw * (1.f - xw)) * sh * (sh - sqmz) / (
          fem::pow2((sh - sqmz)) + fem::pow2(gmmz)) * vint(112) + (
          fem::pow2(vi) + fem::pow2(ai)) / fem::pow2((16.f * xw * (
          1.f - xw))) * sh2 / (fem::pow2((sh - sqmz)) + fem::pow2(
          gmmz)) * vint(114));
        statement_150:;
      }
      //C
    }
    else if (isub == 2) {
      //C...f + fb' -> W+/-.
      pywidt(cmn, 24, fem::sqrt(sh), wdtp, wdte);
      facw = comfac * fem::pow2((aem / xw)) * 1.f / 24 * sh2 / (
        fem::pow2((sh - sqmw)) + fem::pow2(gmmw));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_170;
        }
        ia = fem::iabs(i);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_160;
          }
          ja = fem::iabs(j);
          if (i * j > 0 || fem::mod(ia + ja, 2) == 0) {
            goto statement_160;
          }
          if ((ia <= 10 && ja > 10) || (ia > 10 && ja <= 10)) {
            goto statement_160;
          }
          kchw = (kchg(ia, 1) * fem::isign(1, i) + kchg(ja, 1) * fem::isign(1,
            j)) / 3;
          facf = 1.f;
          if (ia <= 10) {
            facf = vckm((ia + 1) / 2, (ja + 1) / 2) * faca / 3.f;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = facf * facw * (wdte(0, 1) + wdte(0, (5 -
            kchw) / 2) + wdte(0, 4));
          statement_160:;
        }
        statement_170:;
      }
      //C
    }
    else if (isub == 3) {
      //C...f + fb -> H0.
      pywidt(cmn, 25, fem::sqrt(sh), wdtp, wdte);
      fach = comfac * fem::pow2((aem / xw)) * 1.f / 48.f * fem::pow2((
        sh / sqmw)) * sh2 / (fem::pow2((sh - sqmh)) + fem::pow2(
        gmmh)) * (wdte(0, 1) + wdte(0, 2) + wdte(0, 4));
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_180;
        }
        rmq = fem::pow2(pmas(fem::iabs(i), 1)) / sh;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = fach * rmq * fem::sqrt(fem::max(0.f, 1.f - 4.f * rmq));
        statement_180:;
      }
      //C
    }
    else if (isub == 4) {
      //C...gamma + W+/- -> W+/-.
      //C
    }
    else if (isub == 5) {
      //C...Z0 + Z0 -> H0.
      pywidt(cmn, 25, fem::sqrt(sh), wdtp, wdte);
      fach = comfac * 1.f / (128.f * fem::pow2(paru(1)) * 16.f *
        fem::pow3((1.f - xw))) * fem::pow4((aem / xw)) * fem::pow2((
        sh / sqmw)) * sh2 / (fem::pow2((sh - sqmh)) + fem::pow2(
        gmmh)) * (wdte(0, 1) + wdte(0, 2) + wdte(0, 4));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_200;
        }
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_190;
          }
          ei = kchg(fem::iabs(i), 1) / 3.f;
          ai = fem::sign(1.f, ei);
          vi = ai - 4.f * ei * xw;
          ej = kchg(fem::iabs(j), 1) / 3.f;
          aj = fem::sign(1.f, ej);
          vj = aj - 4.f * ej * xw;
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * (fem::pow2(vi) + fem::pow2(ai)) * (
            fem::pow2(vj) + fem::pow2(aj));
          statement_190:;
        }
        statement_200:;
      }
      //C
    }
    else if (isub == 6) {
      //C...Z0 + W+/- -> W+/-.
      //C
    }
    else if (isub == 7) {
      //C...W+ + W- -> Z0.
      //C
    }
    else if (isub == 8) {
      //C...W+ + W- -> H0.
      pywidt(cmn, 25, fem::sqrt(sh), wdtp, wdte);
      fach = comfac * 1.f / (128 * fem::pow2(paru(1))) * fem::pow4((
        aem / xw)) * fem::pow2((sh / sqmw)) * sh2 / (fem::pow2((sh -
        sqmh)) + fem::pow2(gmmh)) * (wdte(0, 1) + wdte(0, 2) + wdte(0,
        4));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_220;
        }
        ei = fem::sign(1.f, fem::ffloat(i)) * kchg(fem::iabs(i), 1);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_210;
          }
          ej = fem::sign(1.f, fem::ffloat(j)) * kchg(fem::iabs(j), 1);
          if (ei * ej > 0.f) {
            goto statement_210;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * vint(180 + i) * vint(180 + j);
          statement_210:;
        }
        statement_220:;
      }
    }
    //C
    //C...B: 2 -> 2, tree diagrams.
    //C
  }
  else if (isub <= 20) {
    if (isub == 11) {
      //C...f + f' -> f + f'.
      facqq1 = comfac * fem::pow2(as) * 4.f / 9.f * (sh2 + uh2) / th2;
      facqqb = comfac * fem::pow2(as) * 4.f / 9.f * ((sh2 + uh2) /
        th2 * faca - mstp(34) * 2.f / 3.f * uh2 / (sh * th));
      facqq2 = comfac * fem::pow2(as) * 4.f / 9.f * ((sh2 + th2) /
        uh2 - mstp(34) * 2.f / 3.f * sh2 / (th * uh));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_240;
        }
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_230;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = facqq1;
          if (i ==  - j) {
            sigh(nchn) = facqqb;
          }
          if (i == j) {
            sigh(nchn) = 0.5f * sigh(nchn);
            nchn++;
            isig(nchn, 1) = i;
            isig(nchn, 2) = j;
            isig(nchn, 3) = 2;
            sigh(nchn) = 0.5f * facqq2;
          }
          statement_230:;
        }
        statement_240:;
      }
      //C
    }
    else if (isub == 12) {
      //C...f + fb -> f' + fb' (q + qb -> q' + qb' only).
      pywidt(cmn, 21, fem::sqrt(sh), wdtp, wdte);
      facqqb = comfac * fem::pow2(as) * 4.f / 9.f * (th2 + uh2) /
        sh2 * (wdte(0, 1) + wdte(0, 2) + wdte(0, 3) + wdte(0, 4));
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_250;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facqqb;
        statement_250:;
      }
      //C
    }
    else if (isub == 13) {
      //C...f + fb -> g + g (q + qb -> g + g only).
      facgg1 = comfac * fem::pow2(as) * 32.f / 27.f * (uh / th - (
        2.f + mstp(34) * 1.f / 4.f) * uh2 / sh2);
      facgg2 = comfac * fem::pow2(as) * 32.f / 27.f * (th / uh - (
        2.f + mstp(34) * 1.f / 4.f) * th2 / sh2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_260;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = 0.5f * facgg1;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 2;
        sigh(nchn) = 0.5f * facgg2;
        statement_260:;
      }
      //C
    }
    else if (isub == 14) {
      //C...f + fb -> g + gamma (q + qb -> g + gamma only).
      facgg = comfac * as * aem * 8.f / 9.f * (th2 + uh2) / (th * uh);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_270;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facgg * fem::pow2(ei);
        statement_270:;
      }
      //C
    }
    else if (isub == 15) {
      //C...f + fb -> g + Z0 (q + qb -> g + Z0 only).
      faczg = comfac * as * aem / (xw * (1.f - xw)) * 1.f / 18.f * (
        th2 + uh2 + 2.f * sqm4 * sh) / (th * uh);
      faczg = faczg * wids(23, 2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_280;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = faczg * (fem::pow2(vi) + fem::pow2(ai));
        statement_280:;
      }
      //C
    }
    else if (isub == 16) {
      //C...f + fb' -> g + W+/- (q + qb' -> g + W+/- only).
      facwg = comfac * as * aem / xw * 2.f / 9.f * (th2 + uh2 + 2.f *
        sqm4 * sh) / (th * uh);
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_300;
        }
        ia = fem::iabs(i);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_290;
          }
          ja = fem::iabs(j);
          if (i * j > 0 || fem::mod(ia + ja, 2) == 0) {
            goto statement_290;
          }
          kchw = (kchg(ia, 1) * fem::isign(1, i) + kchg(ja, 1) * fem::isign(1,
            j)) / 3;
          fckm = 1.f;
          if (mint(43) == 4) {
            fckm = vckm((ia + 1) / 2, (ja + 1) / 2);
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = facwg * fckm * wids(24, (5 - kchw) / 2);
          statement_290:;
        }
        statement_300:;
      }
      //C
    }
    else if (isub == 17) {
      //C...f + fb -> g + H0 (q + qb -> g + H0 only).
      //C
    }
    else if (isub == 18) {
      //C...f + fb -> gamma + gamma.
      facgg = comfac * faca * fem::pow2(aem) * 1.f / 3.f * (th2 +
        uh2) / (th * uh);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_310;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facgg * fem::pow4(ei);
        statement_310:;
      }
      //C
    }
    else if (isub == 19) {
      //C...f + fb -> gamma + Z0.
      facgz = comfac * faca * fem::pow2(aem) / (xw * (1.f - xw)) *
        1.f / 24.f * (th2 + uh2 + 2.f * sqm4 * sh) / (th * uh);
      facgz = facgz * wids(23, 2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_320;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facgz * fem::pow2(ei) * (fem::pow2(vi) + fem::pow2(ai));
        statement_320:;
      }
      //C
    }
    else if (isub == 20) {
      //C...f + fb' -> gamma + W+/-.
      facgw = comfac * faca * fem::pow2(aem) / xw * 1.f / 6.f *
        fem::pow2(((2.f * uh - th) / (3.f * (sh - sqm4)))) * (th2 +
        uh2 + 2.f * sqm4 * sh) / (th * uh);
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_340;
        }
        ia = fem::iabs(i);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_330;
          }
          ja = fem::iabs(j);
          if (i * j > 0 || fem::mod(ia + ja, 2) == 0) {
            goto statement_330;
          }
          kchw = (kchg(ia, 1) * fem::isign(1, i) + kchg(ja, 1) * fem::isign(1,
            j)) / 3;
          fckm = 1.f;
          if (mint(43) == 4) {
            fckm = vckm((ia + 1) / 2, (ja + 1) / 2);
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = facgw * fckm * wids(24, (5 - kchw) / 2);
          statement_330:;
        }
        statement_340:;
      }
    }
    //C
  }
  else if (isub <= 30) {
    if (isub == 21) {
      //C...f + fb -> gamma + H0.
      //C
    }
    else if (isub == 22) {
      //C...f + fb -> Z0 + Z0.
      faczz = comfac * faca * fem::pow2((aem / (xw * (1.f - xw)))) *
        1.f / 768.f * (uh / th + th / uh + 2.f * (sqm3 + sqm4) * sh /
        (th * uh) - sqm3 * sqm4 * (1.f / th2 + 1.f / uh2));
      faczz = faczz * wids(23, 1);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_350;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = faczz * (fem::pow4(vi) + 6.f * fem::pow2(vi) *
          fem::pow2(ai) + fem::pow4(ai));
        statement_350:;
      }
      //C
    }
    else if (isub == 23) {
      //C...f + fb' -> Z0 + W+/-.
      faczw = comfac * faca * fem::pow2((aem / xw)) * 1.f / 6.f;
      faczw = faczw * wids(23, 2);
      thuh = fem::max(th * uh - sqm3 * sqm4, sh * fem::pow2(ckin(3)));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_370;
        }
        ia = fem::iabs(i);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_360;
          }
          ja = fem::iabs(j);
          if (i * j > 0 || fem::mod(ia + ja, 2) == 0) {
            goto statement_360;
          }
          kchw = (kchg(ia, 1) * fem::isign(1, i) + kchg(ja, 1) * fem::isign(1,
            j)) / 3;
          ei = kchg(ia, 1) / 3.f;
          ai = fem::sign(1.f, ei);
          vi = ai - 4.f * ei * xw;
          ej = kchg(ja, 1) / 3.f;
          aj = fem::sign(1.f, ej);
          vj = aj - 4.f * ej * xw;
          if (vi + ai > 0) {
            visav = vi;
            aisav = ai;
            vi = vj;
            ai = aj;
            vj = visav;
            aj = aisav;
          }
          fckm = 1.f;
          if (mint(43) == 4) {
            fckm = vckm((ia + 1) / 2, (ja + 1) / 2);
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = faczw * fckm * (1.f / fem::pow2((sh - sqmw)) *
            ((9.f - 8.f * xw) / 4.f * thuh + (8.f * xw - 6.f) / 4.f *
            sh * (sqm3 + sqm4)) + (thuh - sh * (sqm3 + sqm4)) / (
            2.f * (sh - sqmw)) * ((vj + aj) / th - (vi + ai) / uh) +
            thuh / (16.f * (1.f - xw)) * (fem::pow2((vj + aj)) /
            th2 + fem::pow2((vi + ai)) / uh2) + sh * (sqm3 + sqm4) / (
            8.f * (1.f - xw)) * (vi + ai) * (vj + aj) / (th * uh)) * wids(24,
            (5 - kchw) / 2);
          statement_360:;
        }
        statement_370:;
      }
      //C
    }
    else if (isub == 24) {
      //C...f + fb -> Z0 + H0.
      thuh = fem::max(th * uh - sqm3 * sqm4, sh * fem::pow2(ckin(3)));
      fachz = comfac * faca * fem::pow2((aem / (xw * (1.f - xw)))) *
        1.f / 96.f * (thuh + 2.f * sh * sqmz) / fem::pow2((sh -
        sqmz));
      fachz = fachz * wids(23, 2) * wids(25, 2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_380;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachz * (fem::pow2(vi) + fem::pow2(ai));
        statement_380:;
      }
      //C
    }
    else if (isub == 25) {
      //C...f + fb -> W+ + W-.
      facww = comfac * faca * fem::pow2((aem / xw)) * 1.f / 12.f;
      facww = facww * wids(24, 1);
      thuh = fem::max(th * uh - sqm3 * sqm4, sh * fem::pow2(ckin(3)));
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_390;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        dsigww = thuh / sh2 * (3.f - (sh - 3.f * (sqm3 + sqm4)) / (
          sh - sqmz) * (vi + ai) / (2.f * ai * (1.f - xw)) +
          fem::pow2((sh / (sh - sqmz))) * (1.f - 2.f * (sqm3 +
          sqm4) / sh + 12.f * sqm3 * sqm4 / sh2) * (fem::pow2(vi) +
          fem::pow2(ai)) / (8.f * fem::pow2((1.f - xw)))) - 2.f *
          sqmz / (sh - sqmz) * (vi + ai) / ai + sqmz * sh / fem::pow2(
          (sh - sqmz)) * (1.f - 2.f * (sqm3 + sqm4) / sh) * (
          fem::pow2(vi) + fem::pow2(ai)) / (2.f * (1.f - xw));
        if (kchg(fem::iabs(i), 1) < 0) {
          dsigww += 2.f * (1.f + sqmz / (sh - sqmz) * (vi + ai) / (
            2.f * ai)) * (thuh / (sh * th) - (sqm3 + sqm4) / th) +
            thuh / th2;
        }
        else {
          dsigww += 2.f * (1.f + sqmz / (sh - sqmz) * (vi + ai) / (
            2.f * ai)) * (thuh / (sh * uh) - (sqm3 + sqm4) / uh) +
            thuh / uh2;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facww * dsigww;
        statement_390:;
      }
      //C
    }
    else if (isub == 26) {
      //C...f + fb' -> W+/- + H0.
      thuh = fem::max(th * uh - sqm3 * sqm4, sh * fem::pow2(ckin(3)));
      fachw = comfac * faca * fem::pow2((aem / xw)) * 1.f / 24.f * (
        thuh + 2.f * sh * sqmw) / fem::pow2((sh - sqmw));
      fachw = fachw * wids(25, 2);
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_410;
        }
        ia = fem::iabs(i);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(1, j) == 0) {
            goto statement_400;
          }
          ja = fem::iabs(j);
          if (i * j > 0 || fem::mod(ia + ja, 2) == 0) {
            goto statement_400;
          }
          kchw = (kchg(ia, 1) * fem::isign(1, i) + kchg(ja, 1) * fem::isign(1,
            j)) / 3;
          fckm = 1.f;
          if (mint(43) == 4) {
            fckm = vckm((ia + 1) / 2, (ja + 1) / 2);
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fachw * fckm * wids(24, (5 - kchw) / 2);
          statement_400:;
        }
        statement_410:;
      }
      //C
    }
    else if (isub == 27) {
      //C...f + fb -> H0 + H0.
      //C
    }
    else if (isub == 28) {
      //C...f + g -> f + g (q + g -> q + g only).
      facqg1 = comfac * fem::pow2(as) * 4.f / 9.f * ((2.f + mstp(
        34) * 1.f / 4.f) * uh2 / th2 - uh / sh) * faca;
      facqg2 = comfac * fem::pow2(as) * 4.f / 9.f * ((2.f + mstp(
        34) * 1.f / 4.f) * sh2 / th2 - sh / uh);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0) {
          goto statement_430;
        }
        FEM_DO_SAFE(isde, 1, 2) {
          if (isde == 1 && kfac(1, i) * kfac(2, 21) == 0) {
            goto statement_420;
          }
          if (isde == 2 && kfac(1, 21) * kfac(2, i) == 0) {
            goto statement_420;
          }
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 1;
          sigh(nchn) = facqg1;
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 2;
          sigh(nchn) = facqg2;
          statement_420:;
        }
        statement_430:;
      }
      //C
    }
    else if (isub == 29) {
      //C...f + g -> f + gamma (q + g -> q + gamma only).
      fgq = comfac * faca * as * aem * 1.f / 3.f * (sh2 + uh2) / (-sh * uh);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0) {
          goto statement_450;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        facgq = fgq * fem::pow2(ei);
        FEM_DO_SAFE(isde, 1, 2) {
          if (isde == 1 && kfac(1, i) * kfac(2, 21) == 0) {
            goto statement_440;
          }
          if (isde == 2 && kfac(1, 21) * kfac(2, i) == 0) {
            goto statement_440;
          }
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 1;
          sigh(nchn) = facgq;
          statement_440:;
        }
        statement_450:;
      }
      //C
    }
    else if (isub == 30) {
      //C...f + g -> f + Z0 (q + g -> q + Z0 only).
      fzq = comfac * faca * as * aem / (xw * (1.f - xw)) * 1.f /
        48.f * (sh2 + uh2 + 2.f * sqm4 * th) / (-sh * uh);
      fzq = fzq * wids(23, 2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0) {
          goto statement_470;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        faczq = fzq * (fem::pow2(vi) + fem::pow2(ai));
        FEM_DO_SAFE(isde, 1, 2) {
          if (isde == 1 && kfac(1, i) * kfac(2, 21) == 0) {
            goto statement_460;
          }
          if (isde == 2 && kfac(1, 21) * kfac(2, i) == 0) {
            goto statement_460;
          }
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 1;
          sigh(nchn) = faczq;
          statement_460:;
        }
        statement_470:;
      }
    }
    //C
  }
  else if (isub <= 40) {
    if (isub == 31) {
      //C...f + g -> f' + W+/- (q + g -> q' + W+/- only).
      facwq = comfac * faca * as * aem / xw * 1.f / 12.f * (sh2 +
        uh2 + 2.f * sqm4 * th) / (-sh * uh);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0) {
          goto statement_490;
        }
        ia = fem::iabs(i);
        kchw = fem::isign(1, kchg(ia, 1) * fem::isign(1, i));
        FEM_DO_SAFE(isde, 1, 2) {
          if (isde == 1 && kfac(1, i) * kfac(2, 21) == 0) {
            goto statement_480;
          }
          if (isde == 2 && kfac(1, 21) * kfac(2, i) == 0) {
            goto statement_480;
          }
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 1;
          sigh(nchn) = facwq * vint(180 + i) * wids(24, (5 - kchw) / 2);
          statement_480:;
        }
        statement_490:;
      }
      //C
    }
    else if (isub == 32) {
      //C...f + g -> f + H0 (q + g -> q + H0 only).
      //C
    }
    else if (isub == 33) {
      //C...f + gamma -> f + g (q + gamma -> q + g only).
      //C
    }
    else if (isub == 34) {
      //C...f + gamma -> f + gamma.
      //C
    }
    else if (isub == 35) {
      //C...f + gamma -> f + Z0.
      //C
    }
    else if (isub == 36) {
      //C...f + gamma -> f' + W+/-.
      //C
    }
    else if (isub == 37) {
      //C...f + gamma -> f + H0.
      //C
    }
    else if (isub == 38) {
      //C...f + Z0 -> f + g (q + Z0 -> q + g only).
      //C
    }
    else if (isub == 39) {
      //C...f + Z0 -> f + gamma.
      //C
    }
    else if (isub == 40) {
      //C...f + Z0 -> f + Z0.
    }
    //C
  }
  else if (isub <= 50) {
    if (isub == 41) {
      //C...f + Z0 -> f' + W+/-.
      //C
    }
    else if (isub == 42) {
      //C...f + Z0 -> f + H0.
      //C
    }
    else if (isub == 43) {
      //C...f + W+/- -> f' + g (q + W+/- -> q' + g only).
      //C
    }
    else if (isub == 44) {
      //C...f + W+/- -> f' + gamma.
      //C
    }
    else if (isub == 45) {
      //C...f + W+/- -> f' + Z0.
      //C
    }
    else if (isub == 46) {
      //C...f + W+/- -> f' + W+/-.
      //C
    }
    else if (isub == 47) {
      //C...f + W+/- -> f' + H0.
      //C
    }
    else if (isub == 48) {
      //C...f + H0 -> f + g (q + H0 -> q + g only).
      //C
    }
    else if (isub == 49) {
      //C...f + H0 -> f + gamma.
      //C
    }
    else if (isub == 50) {
      //C...f + H0 -> f + Z0.
    }
    //C
  }
  else if (isub <= 60) {
    if (isub == 51) {
      //C...f + H0 -> f' + W+/-.
      //C
    }
    else if (isub == 52) {
      //C...f + H0 -> f + H0.
      //C
    }
    else if (isub == 53) {
      //C...g + g -> f + fb (g + g -> q + qb only).
      pywidt(cmn, 21, fem::sqrt(sh), wdtp, wdte);
      facqq1 = comfac * fem::pow2(as) * 1.f / 6.f * (uh / th - (2.f +
        mstp(34) * 1.f / 4.f) * uh2 / sh2) * (wdte(0, 1) + wdte(0,
        2) + wdte(0, 3) + wdte(0, 4)) * faca;
      facqq2 = comfac * fem::pow2(as) * 1.f / 6.f * (th / uh - (2.f +
        mstp(34) * 1.f / 4.f) * th2 / sh2) * (wdte(0, 1) + wdte(0,
        2) + wdte(0, 3) + wdte(0, 4)) * faca;
      if (kfac(1, 21) * kfac(2, 21) == 0) {
        goto statement_500;
      }
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 1;
      sigh(nchn) = facqq1;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 2;
      sigh(nchn) = facqq2;
      statement_500:;
      //C
    }
    else if (isub == 54) {
      //C...g + gamma -> f + fb (g + gamma -> q + qb only).
      //C
    }
    else if (isub == 55) {
      //C...g + gamma -> f + fb (g + gamma -> q + qb only).
      //C
    }
    else if (isub == 56) {
      //C...g + gamma -> f + fb (g + gamma -> q + qb only).
      //C
    }
    else if (isub == 57) {
      //C...g + gamma -> f + fb (g + gamma -> q + qb only).
      //C
    }
    else if (isub == 58) {
      //C...gamma + gamma -> f + fb.
      //C
    }
    else if (isub == 59) {
      //C...gamma + Z0 -> f + fb.
      //C
    }
    else if (isub == 60) {
      //C...gamma + W+/- -> f + fb'.
    }
    //C
  }
  else if (isub <= 70) {
    if (isub == 61) {
      //C...gamma + H0 -> f + fb.
      //C
    }
    else if (isub == 62) {
      //C...Z0 + Z0 -> f + fb.
      //C
    }
    else if (isub == 63) {
      //C...Z0 + W+/- -> f + fb'.
      //C
    }
    else if (isub == 64) {
      //C...Z0 + H0 -> f + fb.
      //C
    }
    else if (isub == 65) {
      //C...W+ + W- -> f + fb.
      //C
    }
    else if (isub == 66) {
      //C...W+/- + H0 -> f + fb'.
      //C
    }
    else if (isub == 67) {
      //C...H0 + H0 -> f + fb.
      //C
    }
    else if (isub == 68) {
      //C...g + g -> g + g.
      facgg1 = comfac * fem::pow2(as) * 9.f / 4.f * (sh2 / th2 +
        2.f * sh / th + 3.f + 2.f * th / sh + th2 / sh2) * faca;
      facgg2 = comfac * fem::pow2(as) * 9.f / 4.f * (uh2 / sh2 +
        2.f * uh / sh + 3.f + 2.f * sh / uh + sh2 / uh2) * faca;
      facgg3 = comfac * fem::pow2(as) * 9.f / 4.f * (th2 / uh2 +
        2.f * th / uh + 3 + 2.f * uh / th + uh2 / th2);
      if (kfac(1, 21) * kfac(2, 21) == 0) {
        goto statement_510;
      }
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 1;
      sigh(nchn) = 0.5f * facgg1;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 2;
      sigh(nchn) = 0.5f * facgg2;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 3;
      sigh(nchn) = 0.5f * facgg3;
      statement_510:;
      //C
    }
    else if (isub == 69) {
      //C...gamma + gamma -> W+ + W-.
      //C
    }
    else if (isub == 70) {
      //C...gamma + W+/- -> gamma + W+/-.
    }
    //C
  }
  else if (isub <= 80) {
    if (isub == 71) {
      //C...Z0 + Z0 -> Z0 + Z0.
      be2 = 1.f - 4.f * sqmz / sh;
      th = -0.5f * sh * be2 * (1.f - cth);
      uh = -0.5f * sh * be2 * (1.f + cth);
      shang = 1.f / (1.f - xw) * sqmw / sqmz * fem::pow2((1.f + be2));
      ashre = (sh - sqmh) / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      ashim = -gmmh / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      thang = 1.f / (1.f - xw) * sqmw / sqmz * fem::pow2((be2 - cth));
      athre = (th - sqmh) / (fem::pow2((th - sqmh)) + fem::pow2(gmmh)) * thang;
      athim = -gmmh / (fem::pow2((th - sqmh)) + fem::pow2(gmmh)) * thang;
      uhang = 1.f / (1.f - xw) * sqmw / sqmz * fem::pow2((be2 + cth));
      auhre = (uh - sqmh) / (fem::pow2((uh - sqmh)) + fem::pow2(gmmh)) * uhang;
      auhim = -gmmh / (fem::pow2((uh - sqmh)) + fem::pow2(gmmh)) * uhang;
      fach = 0.5f * comfac * 1.f / (4096.f * fem::pow2(paru(1)) *
        16.f * fem::pow2((1.f - xw))) * fem::pow4((aem / xw)) *
        fem::pow2((sh / sqmw)) * (fem::pow2((ashre + athre +
        auhre)) + fem::pow2((ashim + athim + auhim))) * sqmz / sqmw;
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_530;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        avi = fem::pow2(ai) + fem::pow2(vi);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_520;
          }
          ej = kchg(fem::iabs(j), 1) / 3.f;
          aj = fem::sign(1.f, ej);
          vj = aj - 4.f * ej * xw;
          avj = fem::pow2(aj) + fem::pow2(vj);
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * avi * avj;
          statement_520:;
        }
        statement_530:;
      }
      //C
    }
    else if (isub == 72) {
      //C...Z0 + Z0 -> W+ + W-.
      be2 = fem::sqrt((1.f - 4.f * sqmw / sh) * (1.f - 4.f * sqmz / sh));
      cth2 = fem::pow2(cth);
      th = -0.5f * sh * (1.f - 2.f * (sqmw + sqmz) / sh - be2 * cth);
      uh = -0.5f * sh * (1.f - 2.f * (sqmw + sqmz) / sh + be2 * cth);
      shang = 4.f * fem::sqrt(sqmw / (sqmz * (1.f - xw))) * (1.f -
        2.f * sqmw / sh) * (1.f - 2.f * sqmz / sh);
      ashre = (sh - sqmh) / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      ashim = -gmmh / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      atwre = (1.f - xw) / sqmz * sh / (th - sqmw) * (fem::pow2((cth -
        be2)) * (3.f / 2.f + be2 / 2.f * cth - (sqmw + sqmz) / sh +
        fem::pow2((sqmw - sqmz)) / (sh * sqmw)) + 4.f * ((sqmw +
        sqmz) / sh * (1.f - 3.f * cth2) + 8.f * sqmw * sqmz / sh2 * (
        2.f * cth2 - 1.f) + 4.f * (fem::pow2(sqmw) + fem::pow2(sqmz)) /
        sh2 * cth2 + 2.f * (sqmw + sqmz) / sh * be2 * cth));
      atwim = 0.f;
      auwre = (1.f - xw) / sqmz * sh / (uh - sqmw) * (fem::pow2((cth +
        be2)) * (3.f / 2.f - be2 / 2.f * cth - (sqmw + sqmz) / sh +
        fem::pow2((sqmw - sqmz)) / (sh * sqmw)) + 4.f * ((sqmw +
        sqmz) / sh * (1.f - 3.f * cth2) + 8.f * sqmw * sqmz / sh2 * (
        2.f * cth2 - 1.f) + 4.f * (fem::pow2(sqmw) + fem::pow2(sqmz)) /
        sh2 * cth2 - 2.f * (sqmw + sqmz) / sh * be2 * cth));
      auwim = 0.f;
      a4re = 2.f * (1.f - xw) / sqmz * (3.f - cth2 - 4.f * (sqmw + sqmz) / sh);
      a4im = 0.f;
      fach = comfac * 1.f / (4096.f * fem::pow2(paru(1)) * 16.f *
        fem::pow2((1.f - xw))) * fem::pow4((aem / xw)) * fem::pow2((
        sh / sqmw)) * (fem::pow2((ashre + atwre + auwre + a4re)) +
        fem::pow2((ashim + atwim + auwim + a4im))) * sqmz / sqmw;
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_550;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        avi = fem::pow2(ai) + fem::pow2(vi);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_540;
          }
          ej = kchg(fem::iabs(j), 1) / 3.f;
          aj = fem::sign(1.f, ej);
          vj = aj - 4.f * ej * xw;
          avj = fem::pow2(aj) + fem::pow2(vj);
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * avi * avj;
          statement_540:;
        }
        statement_550:;
      }
      //C
    }
    else if (isub == 73) {
      //C...Z0 + W+/- -> Z0 + W+/-.
      be2 = 1.f - 2.f * (sqmz + sqmw) / sh + fem::pow2(((sqmz - sqmw) / sh));
      ep1 = 1.f + (sqmz - sqmw) / sh;
      ep2 = 1.f - (sqmz - sqmw) / sh;
      th = -0.5f * sh * be2 * (1.f - cth);
      uh = fem::pow2((sqmz - sqmw)) / sh - 0.5f * sh * be2 * (1.f + cth);
      thang = fem::sqrt(sqmw / (sqmz * (1.f - xw))) * (be2 - ep1 *
        cth) * (be2 - ep2 * cth);
      athre = (th - sqmh) / (fem::pow2((th - sqmh)) + fem::pow2(gmmh)) * thang;
      athim = -gmmh / (fem::pow2((th - sqmh)) + fem::pow2(gmmh)) * thang;
      aswre = (1.f - xw) / sqmz * sh / (sh - sqmw) * (-be2 *
        fem::pow4((ep1 + ep2)) * cth + 1.f / 4.f * fem::pow2((be2 +
        ep1 * ep2)) * (fem::pow2((ep1 - ep2)) - 4.f * be2 * cth) +
        2.f * be2 * (be2 + ep1 * ep2) * fem::pow2((ep1 + ep2)) *
        cth - 1.f / 16.f * sh / sqmw * fem::pow2((fem::pow2(ep1) -
        fem::pow2(ep2))) * fem::pow2((be2 + ep1 * ep2)));
      aswim = 0.f;
      auwre = (1.f - xw) / sqmz * sh / (uh - sqmw) * (-be2 * (ep2 +
        ep1 * cth) * (ep1 + ep2 * cth) * (be2 + ep1 * ep2) + be2 * (
        ep2 + ep1 * cth) * (be2 + ep1 * ep2 * cth) * (2.f * ep2 -
        ep2 * cth + ep1) - be2 * fem::pow2((ep2 + ep1 * cth)) * (
        be2 - fem::pow2(ep2) * cth) - 1.f / 8.f * fem::pow2((be2 +
        ep1 * ep2 * cth)) * (fem::pow2((ep1 + ep2)) + 2.f * be2 * (
        1.f - cth)) + 1.f / 32.f * sh / sqmw * fem::pow2((be2 + ep1 *
        ep2 * cth)) * fem::pow2((fem::pow2(ep1) - fem::pow2(ep2))) -
        be2 * (ep1 + ep2 * cth) * (ep2 + ep1 * cth) * (be2 + ep1 *
        ep2) + be2 * (ep1 + ep2 * cth) * (be2 + ep1 * ep2 * cth) * (
        2.f * ep1 - ep1 * cth + ep2) - be2 * fem::pow2((ep1 + ep2 *
        cth)) * (be2 - fem::pow2(ep1) * cth) - 1.f / 8.f * fem::pow2((
        be2 + ep1 * ep2 * cth)) * (fem::pow2((ep1 + ep2)) + 2.f *
        be2 * (1.f - cth)) + 1.f / 32.f * sh / sqmw * fem::pow2((
        be2 + ep1 * ep2 * cth)) * fem::pow2((fem::pow2(ep1) -
        fem::pow2(ep2))));
      auwim = 0.f;
      a4re = (1.f - xw) / sqmz * (fem::pow2(ep1) * fem::pow2(ep2) * (
        fem::pow2(cth) - 1.f) - 2.f * be2 * (fem::pow2(ep1) +
        fem::pow2(ep2) + ep1 * ep2) * cth - 2.f * be2 * ep1 * ep2);
      a4im = 0.f;
      fach = comfac * 1.f / (4096.f * fem::pow2(paru(1)) * 4.f * (
        1.f - xw)) * fem::pow4((aem / xw)) * fem::pow2((sh / sqmw)) *
        (fem::pow2((athre + aswre + auwre + a4re)) + fem::pow2((
        athim + aswim + auwim + a4im))) * fem::sqrt(sqmz / sqmw);
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_570;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        avi = fem::pow2(ai) + fem::pow2(vi);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_560;
          }
          ej = kchg(fem::iabs(j), 1) / 3.f;
          aj = fem::sign(1.f, ej);
          vj = ai - 4.f * ej * xw;
          avj = fem::pow2(aj) + fem::pow2(vj);
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * (avi * vint(180 + j) + vint(180 + i) * avj);
          statement_560:;
        }
        statement_570:;
      }
      //C
    }
    else if (isub == 75) {
      //C...W+ + W- -> gamma + gamma.
      //C
    }
    else if (isub == 76) {
      //C...W+ + W- -> Z0 + Z0.
      be2 = fem::sqrt((1.f - 4.f * sqmw / sh) * (1.f - 4.f * sqmz / sh));
      cth2 = fem::pow2(cth);
      th = -0.5f * sh * (1.f - 2.f * (sqmw + sqmz) / sh - be2 * cth);
      uh = -0.5f * sh * (1.f - 2.f * (sqmw + sqmz) / sh + be2 * cth);
      shang = 4.f * fem::sqrt(sqmw / (sqmz * (1.f - xw))) * (1.f -
        2.f * sqmw / sh) * (1.f - 2.f * sqmz / sh);
      ashre = (sh - sqmh) / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      ashim = -gmmh / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      atwre = (1.f - xw) / sqmz * sh / (th - sqmw) * (fem::pow2((cth -
        be2)) * (3.f / 2.f + be2 / 2.f * cth - (sqmw + sqmz) / sh +
        fem::pow2((sqmw - sqmz)) / (sh * sqmw)) + 4.f * ((sqmw +
        sqmz) / sh * (1.f - 3.f * cth2) + 8.f * sqmw * sqmz / sh2 * (
        2.f * cth2 - 1.f) + 4.f * (fem::pow2(sqmw) + fem::pow2(sqmz)) /
        sh2 * cth2 + 2.f * (sqmw + sqmz) / sh * be2 * cth));
      atwim = 0.f;
      auwre = (1.f - xw) / sqmz * sh / (uh - sqmw) * (fem::pow2((cth +
        be2)) * (3.f / 2.f - be2 / 2.f * cth - (sqmw + sqmz) / sh +
        fem::pow2((sqmw - sqmz)) / (sh * sqmw)) + 4.f * ((sqmw +
        sqmz) / sh * (1.f - 3.f * cth2) + 8.f * sqmw * sqmz / sh2 * (
        2.f * cth2 - 1.f) + 4.f * (fem::pow2(sqmw) + fem::pow2(sqmz)) /
        sh2 * cth2 - 2.f * (sqmw + sqmz) / sh * be2 * cth));
      auwim = 0.f;
      a4re = 2.f * (1.f - xw) / sqmz * (3.f - cth2 - 4.f * (sqmw + sqmz) / sh);
      a4im = 0.f;
      fach = 0.5f * comfac * 1.f / (4096.f * fem::pow2(paru(1))) *
        fem::pow4((aem / xw)) * fem::pow2((sh / sqmw)) * (fem::pow2((
        ashre + atwre + auwre + a4re)) + fem::pow2((ashim + atwim +
        auwim + a4im)));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_590;
        }
        ei = fem::sign(1.f, fem::ffloat(i)) * kchg(fem::iabs(i), 1);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_580;
          }
          ej = fem::sign(1.f, fem::ffloat(j)) * kchg(fem::iabs(j), 1);
          if (ei * ej > 0.f) {
            goto statement_580;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * vint(180 + i) * vint(180 + j);
          statement_580:;
        }
        statement_590:;
      }
      //C
    }
    else if (isub == 77) {
      //C...W+/- + W+/- -> W+/- + W+/-.
      be2 = 1.f - 4.f * sqmw / sh;
      be4 = fem::pow2(be2);
      cth2 = fem::pow2(cth);
      cth3 = fem::pow3(cth);
      th = -0.5f * sh * be2 * (1.f - cth);
      uh = -0.5f * sh * be2 * (1.f + cth);
      shang = fem::pow2((1.f + be2));
      ashre = (sh - sqmh) / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      ashim = -gmmh / (fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * shang;
      thang = fem::pow2((be2 - cth));
      athre = (th - sqmh) / (fem::pow2((th - sqmh)) + fem::pow2(gmmh)) * thang;
      athim = -gmmh / (fem::pow2((th - sqmh)) + fem::pow2(gmmh)) * thang;
      sgzang = 1.f / sqmw * be2 * fem::pow2((3.f - be2)) * cth;
      asgre = xw * sgzang;
      asgim = 0.f;
      aszre = (1.f - xw) * sh / (sh - sqmz) * sgzang;
      aszim = 0.f;
      tgzang = 1.f / sqmw * (be2 * (4.f - 2.f * be2 + be4) + be2 * (
        4.f - 10.f * be2 + be4) * cth + (2.f - 11.f * be2 + 10.f *
        be4) * cth2 + be2 * cth3);
      atgre = 0.5f * xw * sh / th * tgzang;
      atgim = 0.f;
      atzre = 0.5f * (1.f - xw) * sh / (th - sqmz) * tgzang;
      atzim = 0.f;
      a4re = 1.f / sqmw * (1.f + 2.f * be2 - 6.f * be2 * cth - cth2);
      a4im = 0.f;
      fach = comfac * 1.f / (4096.f * fem::pow2(paru(1))) * fem::pow4(
        (aem / xw)) * fem::pow2((sh / sqmw)) * (fem::pow2((ashre +
        athre + asgre + aszre + atgre + atzre + a4re)) + fem::pow2((
        ashim + athim + asgim + aszim + atgim + atzim + a4im)));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_610;
        }
        ei = fem::sign(1.f, fem::ffloat(i)) * kchg(fem::iabs(i), 1);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_600;
          }
          ej = fem::sign(1.f, fem::ffloat(j)) * kchg(fem::iabs(j), 1);
          if (ei * ej > 0.f) {
            goto statement_600;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = fach * vint(180 + i) * vint(180 + j);
          statement_600:;
        }
        statement_610:;
      }
      //C
    }
    else if (isub == 78) {
      //C...W+/- + H0 -> W+/- + H0.
      //C
    }
    else if (isub == 79) {
      //C...H0 + H0 -> H0 + H0.
      //C
    }
    //C
    //C...C: 2 -> 2, tree diagrams with masses.
    //C
  }
  else if (isub <= 90) {
    if (isub == 81) {
      //C...q + qb -> Q + QB.
      facqqb = comfac * fem::pow2(as) * 4.f / 9.f * ((fem::pow2((th -
        sqm3)) + fem::pow2((uh - sqm3))) / sh2 + 2.f * sqm3 / sh);
      if (mstp(35) >= 1) {
        if (mstp(35) == 1) {
          alssg = parp(35);
        }
        else {
          mst115 = mstu(115);
          mstu(115) = mstp(36);
          q2bn = fem::sqrt(sqm3 * (fem::pow2((fem::sqrt(sh) - 2.f *
            fem::sqrt(sqm3))) + fem::pow2(parp(36))));
          alssg = ulalps(cmn, q2bn);
          mstu(115) = mst115;
        }
        xrepu = paru(1) * alssg / (6.f * fem::sqrt(fem::max(1e-20f,
          1.f - 4.f * sqm3 / sh)));
        frepu = xrepu / (fem::exp(fem::min(100.f, xrepu)) - 1.f);
        pari(81) = frepu;
        facqqb = facqqb * frepu;
      }
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_620;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facqqb;
        statement_620:;
      }
      //C
    }
    else if (isub == 82) {
      //C...g + g -> Q + QB.
      facqq1 = comfac * faca * fem::pow2(as) * 1.f / 6.f * ((uh -
        sqm3) / (th - sqm3) - 2.f * fem::pow2((uh - sqm3)) / sh2 +
        4.f * sqm3 / sh * (th * uh - fem::pow2(sqm3)) / fem::pow2((
        th - sqm3)));
      facqq2 = comfac * faca * fem::pow2(as) * 1.f / 6.f * ((th -
        sqm3) / (uh - sqm3) - 2.f * fem::pow2((th - sqm3)) / sh2 +
        4.f * sqm3 / sh * (th * uh - fem::pow2(sqm3)) / fem::pow2((
        uh - sqm3)));
      if (mstp(35) >= 1) {
        if (mstp(35) == 1) {
          alssg = parp(35);
        }
        else {
          mst115 = mstu(115);
          mstu(115) = mstp(36);
          q2bn = fem::sqrt(sqm3 * (fem::pow2((fem::sqrt(sh) - 2.f *
            fem::sqrt(sqm3))) + fem::pow2(parp(36))));
          alssg = ulalps(cmn, q2bn);
          mstu(115) = mst115;
        }
        xattr = 4.f * paru(1) * alssg / (3.f * fem::sqrt(fem::max(1e-20f,
          1.f - 4.f * sqm3 / sh)));
        fattr = xattr / (1.f - fem::exp(-fem::min(100.f, xattr)));
        xrepu = paru(1) * alssg / (6.f * fem::sqrt(fem::max(1e-20f,
          1.f - 4.f * sqm3 / sh)));
        frepu = xrepu / (fem::exp(fem::min(100.f, xrepu)) - 1.f);
        fatre = (2.f * fattr + 5.f * frepu) / 7.f;
        pari(81) = fatre;
        facqq1 = facqq1 * fatre;
        facqq2 = facqq2 * fatre;
      }
      if (kfac(1, 21) * kfac(2, 21) == 0) {
        goto statement_630;
      }
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 1;
      sigh(nchn) = facqq1;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 2;
      sigh(nchn) = facqq2;
      statement_630:;
      //C
    }
    //C
    //C...D: Mimimum bias processes.
    //C
  }
  else if (isub <= 100) {
    if (isub == 91) {
      //C...Elastic scattering.
      sigs = xsec(isub, 1);
      //C
    }
    else if (isub == 92) {
      //C...Single diffractive scattering.
      sigs = xsec(isub, 1);
      //C
    }
    else if (isub == 93) {
      //C...Double diffractive scattering.
      sigs = xsec(isub, 1);
      //C
    }
    else if (isub == 94) {
      //C...Central diffractive scattering.
      sigs = xsec(isub, 1);
      //C
    }
    else if (isub == 95) {
      //C...Low-pT scattering.
      sigs = xsec(isub, 1);
      //C
    }
    else if (isub == 96) {
      //C...Multiple interactions: sum of QCD processes.
      pywidt(cmn, 21, fem::sqrt(sh), wdtp, wdte);
      //C
      //C...q + q' -> q + q'.
      facqq1 = comfac * fem::pow2(as) * 4.f / 9.f * (sh2 + uh2) / th2;
      facqqb = comfac * fem::pow2(as) * 4.f / 9.f * ((sh2 + uh2) /
        th2 * faca - mstp(34) * 2.f / 3.f * uh2 / (sh * th));
      facqq2 = comfac * fem::pow2(as) * 4.f / 9.f * ((sh2 + th2) /
        uh2 - mstp(34) * 2.f / 3.f * sh2 / (th * uh));
      FEM_DO_SAFE(i, -3, 3) {
        if (i == 0) {
          goto statement_650;
        }
        FEM_DO_SAFE(j, -3, 3) {
          if (j == 0) {
            goto statement_640;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 111;
          sigh(nchn) = facqq1;
          if (i ==  - j) {
            sigh(nchn) = facqqb;
          }
          if (i == j) {
            sigh(nchn) = 0.5f * sigh(nchn);
            nchn++;
            isig(nchn, 1) = i;
            isig(nchn, 2) = j;
            isig(nchn, 3) = 112;
            sigh(nchn) = 0.5f * facqq2;
          }
          statement_640:;
        }
        statement_650:;
      }
      //C
      //C...q + qb -> q' + qb' or g + g.
      facqqb = comfac * fem::pow2(as) * 4.f / 9.f * (th2 + uh2) /
        sh2 * (wdte(0, 1) + wdte(0, 2) + wdte(0, 3) + wdte(0, 4));
      facgg1 = comfac * fem::pow2(as) * 32.f / 27.f * (uh / th - (
        2.f + mstp(34) * 1.f / 4.f) * uh2 / sh2);
      facgg2 = comfac * fem::pow2(as) * 32.f / 27.f * (th / uh - (
        2.f + mstp(34) * 1.f / 4.f) * th2 / sh2);
      FEM_DO_SAFE(i, -3, 3) {
        if (i == 0) {
          goto statement_660;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 121;
        sigh(nchn) = facqqb;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 131;
        sigh(nchn) = 0.5f * facgg1;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 132;
        sigh(nchn) = 0.5f * facgg2;
        statement_660:;
      }
      //C
      //C...q + g -> q + g.
      facqg1 = comfac * fem::pow2(as) * 4.f / 9.f * ((2.f + mstp(
        34) * 1.f / 4.f) * uh2 / th2 - uh / sh) * faca;
      facqg2 = comfac * fem::pow2(as) * 4.f / 9.f * ((2.f + mstp(
        34) * 1.f / 4.f) * sh2 / th2 - sh / uh);
      FEM_DO_SAFE(i, -3, 3) {
        if (i == 0) {
          goto statement_680;
        }
        FEM_DO_SAFE(isde, 1, 2) {
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 281;
          sigh(nchn) = facqg1;
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 282;
          sigh(nchn) = facqg2;
        }
        statement_680:;
      }
      //C
      //C...g + g -> q + qb or g + g.
      facqq1 = comfac * fem::pow2(as) * 1.f / 6.f * (uh / th - (2.f +
        mstp(34) * 1.f / 4.f) * uh2 / sh2) * (wdte(0, 1) + wdte(0,
        2) + wdte(0, 3) + wdte(0, 4)) * faca;
      facqq2 = comfac * fem::pow2(as) * 1.f / 6.f * (th / uh - (2.f +
        mstp(34) * 1.f / 4.f) * th2 / sh2) * (wdte(0, 1) + wdte(0,
        2) + wdte(0, 3) + wdte(0, 4)) * faca;
      facgg1 = comfac * fem::pow2(as) * 9.f / 4.f * (sh2 / th2 +
        2.f * sh / th + 3.f + 2.f * th / sh + th2 / sh2) * faca;
      facgg2 = comfac * fem::pow2(as) * 9.f / 4.f * (uh2 / sh2 +
        2.f * uh / sh + 3.f + 2.f * sh / uh + sh2 / uh2) * faca;
      facgg3 = comfac * fem::pow2(as) * 9.f / 4.f * (th2 / uh2 +
        2.f * th / uh + 3 + 2.f * uh / th + uh2 / th2);
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 531;
      sigh(nchn) = facqq1;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 532;
      sigh(nchn) = facqq2;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 681;
      sigh(nchn) = 0.5f * facgg1;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 682;
      sigh(nchn) = 0.5f * facgg2;
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 683;
      sigh(nchn) = 0.5f * facgg3;
    }
    //C
    //C...E: 2 -> 1, loop diagrams.
    //C
  }
  else if (isub <= 110) {
    if (isub == 101) {
      //C...g + g -> gamma*/Z0.
      //C
    }
    else if (isub == 102) {
      //C...g + g -> H0.
      pywidt(cmn, 25, fem::sqrt(sh), wdtp, wdte);
      etare = 0.f;
      etaim = 0.f;
      FEM_DO_SAFE(i, 1, 2 * mstp(1)) {
        eps = 4.f * fem::pow2(pmas(i, 1)) / sh;
        if (eps <= 1.f) {
          if (eps > 1.e-4f) {
            root = fem::sqrt(1.f - eps);
            rln = fem::log((1.f + root) / (1.f - root));
          }
          else {
            rln = fem::log(4.f / eps - 2.f);
          }
          phire = 0.25f * (fem::pow2(rln) - fem::pow2(paru(1)));
          phiim = 0.5f * paru(1) * rln;
        }
        else {
          phire = -fem::pow2((fem::asin(1.f / fem::sqrt(eps))));
          phiim = 0.f;
        }
        etare += 0.5f * eps * (1.f + (eps - 1.f) * phire);
        etaim += 0.5f * eps * (eps - 1.f) * phiim;
      }
      eta2 = fem::pow2(etare) + fem::pow2(etaim);
      fach = comfac * faca * fem::pow2((as / paru(1) * aem / xw)) *
        1.f / 512.f * fem::pow2((sh / sqmw)) * eta2 * sh2 / (
        fem::pow2((sh - sqmh)) + fem::pow2(gmmh)) * (wdte(0, 1) + wdte(0,
        2) + wdte(0, 4));
      if (kfac(1, 21) * kfac(2, 21) == 0) {
        goto statement_700;
      }
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 1;
      sigh(nchn) = fach;
      statement_700:;
      //C
    }
    //C
    //C...F: 2 -> 2, box diagrams.
    //C
  }
  else if (isub <= 120) {
    if (isub == 111) {
      //C...f + fb -> g + H0 (q + qb -> g + H0 only).
      a5stur = 0.f;
      a5stui = 0.f;
      FEM_DO_SAFE(i, 1, 2 * mstp(1)) {
        sqmq = fem::pow2(pmas(i, 1));
        epss = 4.f * sqmq / sh;
        epsh = 4.f * sqmq / sqmh;
        a5stur += sqmq / sqmh * (4.f + 4.f * sh / (th + uh) * (pyw1au(cmn,
          epss, 1) - pyw1au(cmn, epsh, 1)) + (1.f - 4.f * sqmq / (
          th + uh)) * (pyw2au(cmn, epss, 1) - pyw2au(cmn, epsh, 1)));
        a5stui += sqmq / sqmh * (4.f * sh / (th + uh) * (pyw1au(cmn,
          epss, 2) - pyw1au(cmn, epsh, 2)) + (1.f - 4.f * sqmq / (
          th + uh)) * (pyw2au(cmn, epss, 2) - pyw2au(cmn, epsh, 2)));
      }
      facgh = comfac * faca / (144.f * fem::pow2(paru(1))) * aem /
        xw * fem::pow3(as) * sqmh / sqmw * sqmh / sh * (fem::pow2(
        uh) + fem::pow2(th)) / fem::pow2((uh + th)) * (fem::pow2(
        a5stur) + fem::pow2(a5stui));
      facgh = facgh * wids(25, 2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_720;
        }
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = facgh;
        statement_720:;
      }
      //C
    }
    else if (isub == 112) {
      //C...f + g -> f + H0 (q + g -> q + H0 only).
      a5tsur = 0.f;
      a5tsui = 0.f;
      FEM_DO_SAFE(i, 1, 2 * mstp(1)) {
        sqmq = fem::pow2(pmas(i, 1));
        epst = 4.f * sqmq / th;
        epsh = 4.f * sqmq / sqmh;
        a5tsur += sqmq / sqmh * (4.f + 4.f * th / (sh + uh) * (pyw1au(cmn,
          epst, 1) - pyw1au(cmn, epsh, 1)) + (1.f - 4.f * sqmq / (
          sh + uh)) * (pyw2au(cmn, epst, 1) - pyw2au(cmn, epsh, 1)));
        a5tsui += sqmq / sqmh * (4.f * th / (sh + uh) * (pyw1au(cmn,
          epst, 2) - pyw1au(cmn, epsh, 2)) + (1.f - 4.f * sqmq / (
          sh + uh)) * (pyw2au(cmn, epst, 2) - pyw2au(cmn, epsh, 2)));
      }
      facqh = comfac * faca / (384.f * fem::pow2(paru(1))) * aem /
        xw * fem::pow3(as) * sqmh / sqmw * sqmh / (-th) * (fem::pow2(
        uh) + fem::pow2(sh)) / fem::pow2((uh + sh)) * (fem::pow2(
        a5tsur) + fem::pow2(a5tsui));
      facqh = facqh * wids(25, 2);
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0) {
          goto statement_750;
        }
        FEM_DO_SAFE(isde, 1, 2) {
          if (isde == 1 && kfac(1, i) * kfac(2, 21) == 0) {
            goto statement_740;
          }
          if (isde == 2 && kfac(1, 21) * kfac(2, i) == 0) {
            goto statement_740;
          }
          nchn++;
          isig(nchn, isde) = i;
          isig(nchn, 3 - isde) = 21;
          isig(nchn, 3) = 1;
          sigh(nchn) = facqh;
          statement_740:;
        }
        statement_750:;
      }
      //C
    }
    else if (isub == 113) {
      //C...g + g -> g + H0.
      a2stur = 0.f;
      a2stui = 0.f;
      a2ustr = 0.f;
      a2usti = 0.f;
      a2tusr = 0.f;
      a2tusi = 0.f;
      a4stur = 0.f;
      a4stui = 0.f;
      FEM_DO_SAFE(i, 6, 2 * mstp(1)) {
        //C'''Only t-quarks yet included
        sqmq = fem::pow2(pmas(i, 1));
        epss = 4.f * sqmq / sh;
        epst = 4.f * sqmq / th;
        epsu = 4.f * sqmq / uh;
        epsh = 4.f * sqmq / sqmh;
        if (epsh < 1.e-6f) {
          goto statement_760;
        }
        bestu = 0.5f * (1.f + fem::sqrt(1.f + epss * th / uh));
        beust = 0.5f * (1.f + fem::sqrt(1.f + epsu * sh / th));
        betus = 0.5f * (1.f + fem::sqrt(1.f + epst * uh / sh));
        beuts = bestu;
        betsu = beust;
        besut = betus;
        w3stur = pyi3au(cmn, bestu, epsh, 1) - pyi3au(cmn, bestu,
          epss, 1) - pyi3au(cmn, bestu, epsu, 1);
        w3stui = pyi3au(cmn, bestu, epsh, 2) - pyi3au(cmn, bestu,
          epss, 2) - pyi3au(cmn, bestu, epsu, 2);
        w3sutr = pyi3au(cmn, besut, epsh, 1) - pyi3au(cmn, besut,
          epss, 1) - pyi3au(cmn, besut, epst, 1);
        w3suti = pyi3au(cmn, besut, epsh, 2) - pyi3au(cmn, besut,
          epss, 2) - pyi3au(cmn, besut, epst, 2);
        w3tsur = pyi3au(cmn, betsu, epsh, 1) - pyi3au(cmn, betsu,
          epst, 1) - pyi3au(cmn, betsu, epsu, 1);
        w3tsui = pyi3au(cmn, betsu, epsh, 2) - pyi3au(cmn, betsu,
          epst, 2) - pyi3au(cmn, betsu, epsu, 2);
        w3tusr = pyi3au(cmn, betus, epsh, 1) - pyi3au(cmn, betus,
          epst, 1) - pyi3au(cmn, betus, epss, 1);
        w3tusi = pyi3au(cmn, betus, epsh, 2) - pyi3au(cmn, betus,
          epst, 2) - pyi3au(cmn, betus, epss, 2);
        w3ustr = pyi3au(cmn, beust, epsh, 1) - pyi3au(cmn, beust,
          epsu, 1) - pyi3au(cmn, beust, epst, 1);
        w3usti = pyi3au(cmn, beust, epsh, 2) - pyi3au(cmn, beust,
          epsu, 2) - pyi3au(cmn, beust, epst, 2);
        w3utsr = pyi3au(cmn, beuts, epsh, 1) - pyi3au(cmn, beuts,
          epsu, 1) - pyi3au(cmn, beuts, epss, 1);
        w3utsi = pyi3au(cmn, beuts, epsh, 2) - pyi3au(cmn, beuts,
          epsu, 2) - pyi3au(cmn, beuts, epss, 2);
        b2stur = sqmq / fem::pow2(sqmh) * (sh * (uh - sh) / (sh +
          uh) + 2.f * th * uh * (uh + 2.f * sh) / fem::pow2((sh +
          uh)) * (pyw1au(cmn, epst, 1) - pyw1au(cmn, epsh, 1)) + (
          sqmq - sh / 4.f) * (0.5f * pyw2au(cmn, epss, 1) + 0.5f * pyw2au(cmn,
          epsh, 1) - pyw2au(cmn, epst, 1) + w3stur) + fem::pow2(sh) *
          (2.f * sqmq / fem::pow2((sh + uh)) - 0.5f / (sh + uh)) * (pyw2au(cmn,
          epst, 1) - pyw2au(cmn, epsh, 1)) + 0.5f * th * uh / sh * (pyw2au(cmn,
          epsh, 1) - 2.f * pyw2au(cmn, epst, 1)) + 0.125f * (sh -
          12.f * sqmq - 4.f * th * uh / sh) * w3tsur);
        b2stui = sqmq / fem::pow2(sqmh) * (2.f * th * uh * (uh +
          2.f * sh) / fem::pow2((sh + uh)) * (pyw1au(cmn, epst, 2) -
          pyw1au(cmn, epsh, 2)) + (sqmq - sh / 4.f) * (0.5f * pyw2au(cmn,
          epss, 2) + 0.5f * pyw2au(cmn, epsh, 2) - pyw2au(cmn, epst,
          2) + w3stui) + fem::pow2(sh) * (2.f * sqmq / fem::pow2((
          sh + uh)) - 0.5f / (sh + uh)) * (pyw2au(cmn, epst, 2) - pyw2au(cmn,
          epsh, 2)) + 0.5f * th * uh / sh * (pyw2au(cmn, epsh, 2) -
          2.f * pyw2au(cmn, epst, 2)) + 0.125f * (sh - 12.f * sqmq -
          4.f * th * uh / sh) * w3tsui);
        b2sutr = sqmq / fem::pow2(sqmh) * (sh * (th - sh) / (sh +
          th) + 2.f * uh * th * (th + 2.f * sh) / fem::pow2((sh +
          th)) * (pyw1au(cmn, epsu, 1) - pyw1au(cmn, epsh, 1)) + (
          sqmq - sh / 4.f) * (0.5f * pyw2au(cmn, epss, 1) + 0.5f * pyw2au(cmn,
          epsh, 1) - pyw2au(cmn, epsu, 1) + w3sutr) + fem::pow2(sh) *
          (2.f * sqmq / fem::pow2((sh + th)) - 0.5f / (sh + th)) * (pyw2au(cmn,
          epsu, 1) - pyw2au(cmn, epsh, 1)) + 0.5f * uh * th / sh * (pyw2au(cmn,
          epsh, 1) - 2.f * pyw2au(cmn, epsu, 1)) + 0.125f * (sh -
          12.f * sqmq - 4.f * uh * th / sh) * w3ustr);
        b2suti = sqmq / fem::pow2(sqmh) * (2.f * uh * th * (th +
          2.f * sh) / fem::pow2((sh + th)) * (pyw1au(cmn, epsu, 2) -
          pyw1au(cmn, epsh, 2)) + (sqmq - sh / 4.f) * (0.5f * pyw2au(cmn,
          epss, 2) + 0.5f * pyw2au(cmn, epsh, 2) - pyw2au(cmn, epsu,
          2) + w3suti) + fem::pow2(sh) * (2.f * sqmq / fem::pow2((
          sh + th)) - 0.5f / (sh + th)) * (pyw2au(cmn, epsu, 2) - pyw2au(cmn,
          epsh, 2)) + 0.5f * uh * th / sh * (pyw2au(cmn, epsh, 2) -
          2.f * pyw2au(cmn, epsu, 2)) + 0.125f * (sh - 12.f * sqmq -
          4.f * uh * th / sh) * w3usti);
        b2tsur = sqmq / fem::pow2(sqmh) * (th * (uh - th) / (th +
          uh) + 2.f * sh * uh * (uh + 2.f * th) / fem::pow2((th +
          uh)) * (pyw1au(cmn, epss, 1) - pyw1au(cmn, epsh, 1)) + (
          sqmq - th / 4.f) * (0.5f * pyw2au(cmn, epst, 1) + 0.5f * pyw2au(cmn,
          epsh, 1) - pyw2au(cmn, epss, 1) + w3tsur) + fem::pow2(th) *
          (2.f * sqmq / fem::pow2((th + uh)) - 0.5f / (th + uh)) * (pyw2au(cmn,
          epss, 1) - pyw2au(cmn, epsh, 1)) + 0.5f * sh * uh / th * (pyw2au(cmn,
          epsh, 1) - 2.f * pyw2au(cmn, epss, 1)) + 0.125f * (th -
          12.f * sqmq - 4.f * sh * uh / th) * w3stur);
        b2tsui = sqmq / fem::pow2(sqmh) * (2.f * sh * uh * (uh +
          2.f * th) / fem::pow2((th + uh)) * (pyw1au(cmn, epss, 2) -
          pyw1au(cmn, epsh, 2)) + (sqmq - th / 4.f) * (0.5f * pyw2au(cmn,
          epst, 2) + 0.5f * pyw2au(cmn, epsh, 2) - pyw2au(cmn, epss,
          2) + w3tsui) + fem::pow2(th) * (2.f * sqmq / fem::pow2((
          th + uh)) - 0.5f / (th + uh)) * (pyw2au(cmn, epss, 2) - pyw2au(cmn,
          epsh, 2)) + 0.5f * sh * uh / th * (pyw2au(cmn, epsh, 2) -
          2.f * pyw2au(cmn, epss, 2)) + 0.125f * (th - 12.f * sqmq -
          4.f * sh * uh / th) * w3stui);
        b2tusr = sqmq / fem::pow2(sqmh) * (th * (sh - th) / (th +
          sh) + 2.f * uh * sh * (sh + 2.f * th) / fem::pow2((th +
          sh)) * (pyw1au(cmn, epsu, 1) - pyw1au(cmn, epsh, 1)) + (
          sqmq - th / 4.f) * (0.5f * pyw2au(cmn, epst, 1) + 0.5f * pyw2au(cmn,
          epsh, 1) - pyw2au(cmn, epsu, 1) + w3tusr) + fem::pow2(th) *
          (2.f * sqmq / fem::pow2((th + sh)) - 0.5f / (th + sh)) * (pyw2au(cmn,
          epsu, 1) - pyw2au(cmn, epsh, 1)) + 0.5f * uh * sh / th * (pyw2au(cmn,
          epsh, 1) - 2.f * pyw2au(cmn, epsu, 1)) + 0.125f * (th -
          12.f * sqmq - 4.f * uh * sh / th) * w3utsr);
        b2tusi = sqmq / fem::pow2(sqmh) * (2.f * uh * sh * (sh +
          2.f * th) / fem::pow2((th + sh)) * (pyw1au(cmn, epsu, 2) -
          pyw1au(cmn, epsh, 2)) + (sqmq - th / 4.f) * (0.5f * pyw2au(cmn,
          epst, 2) + 0.5f * pyw2au(cmn, epsh, 2) - pyw2au(cmn, epsu,
          2) + w3tusi) + fem::pow2(th) * (2.f * sqmq / fem::pow2((
          th + sh)) - 0.5f / (th + sh)) * (pyw2au(cmn, epsu, 2) - pyw2au(cmn,
          epsh, 2)) + 0.5f * uh * sh / th * (pyw2au(cmn, epsh, 2) -
          2.f * pyw2au(cmn, epsu, 2)) + 0.125f * (th - 12.f * sqmq -
          4.f * uh * sh / th) * w3utsi);
        b2ustr = sqmq / fem::pow2(sqmh) * (uh * (th - uh) / (uh +
          th) + 2.f * sh * th * (th + 2.f * uh) / fem::pow2((uh +
          th)) * (pyw1au(cmn, epss, 1) - pyw1au(cmn, epsh, 1)) + (
          sqmq - uh / 4.f) * (0.5f * pyw2au(cmn, epsu, 1) + 0.5f * pyw2au(cmn,
          epsh, 1) - pyw2au(cmn, epss, 1) + w3ustr) + fem::pow2(uh) *
          (2.f * sqmq / fem::pow2((uh + th)) - 0.5f / (uh + th)) * (pyw2au(cmn,
          epss, 1) - pyw2au(cmn, epsh, 1)) + 0.5f * sh * th / uh * (pyw2au(cmn,
          epsh, 1) - 2.f * pyw2au(cmn, epss, 1)) + 0.125f * (uh -
          12.f * sqmq - 4.f * sh * th / uh) * w3sutr);
        b2usti = sqmq / fem::pow2(sqmh) * (2.f * sh * th * (th +
          2.f * uh) / fem::pow2((uh + th)) * (pyw1au(cmn, epss, 2) -
          pyw1au(cmn, epsh, 2)) + (sqmq - uh / 4.f) * (0.5f * pyw2au(cmn,
          epsu, 2) + 0.5f * pyw2au(cmn, epsh, 2) - pyw2au(cmn, epss,
          2) + w3usti) + fem::pow2(uh) * (2.f * sqmq / fem::pow2((
          uh + th)) - 0.5f / (uh + th)) * (pyw2au(cmn, epss, 2) - pyw2au(cmn,
          epsh, 2)) + 0.5f * sh * th / uh * (pyw2au(cmn, epsh, 2) -
          2.f * pyw2au(cmn, epss, 2)) + 0.125f * (uh - 12.f * sqmq -
          4.f * sh * th / uh) * w3suti);
        b2utsr = sqmq / fem::pow2(sqmh) * (uh * (sh - uh) / (uh +
          sh) + 2.f * th * sh * (sh + 2.f * uh) / fem::pow2((uh +
          sh)) * (pyw1au(cmn, epst, 1) - pyw1au(cmn, epsh, 1)) + (
          sqmq - uh / 4.f) * (0.5f * pyw2au(cmn, epsu, 1) + 0.5f * pyw2au(cmn,
          epsh, 1) - pyw2au(cmn, epst, 1) + w3utsr) + fem::pow2(uh) *
          (2.f * sqmq / fem::pow2((uh + sh)) - 0.5f / (uh + sh)) * (pyw2au(cmn,
          epst, 1) - pyw2au(cmn, epsh, 1)) + 0.5f * th * sh / uh * (pyw2au(cmn,
          epsh, 1) - 2.f * pyw2au(cmn, epst, 1)) + 0.125f * (uh -
          12.f * sqmq - 4.f * th * sh / uh) * w3tusr);
        b2utsi = sqmq / fem::pow2(sqmh) * (2.f * th * sh * (sh +
          2.f * uh) / fem::pow2((uh + sh)) * (pyw1au(cmn, epst, 2) -
          pyw1au(cmn, epsh, 2)) + (sqmq - uh / 4.f) * (0.5f * pyw2au(cmn,
          epsu, 2) + 0.5f * pyw2au(cmn, epsh, 2) - pyw2au(cmn, epst,
          2) + w3utsi) + fem::pow2(uh) * (2.f * sqmq / fem::pow2((
          uh + sh)) - 0.5f / (uh + sh)) * (pyw2au(cmn, epst, 2) - pyw2au(cmn,
          epsh, 2)) + 0.5f * th * sh / uh * (pyw2au(cmn, epsh, 2) -
          2.f * pyw2au(cmn, epst, 2)) + 0.125f * (uh - 12.f * sqmq -
          4.f * th * sh / uh) * w3tusi);
        b4stur = sqmq / sqmh * (-2.f / 3.f + (sqmq / sqmh - 1.f /
          4.f) * (pyw2au(cmn, epss, 1) - pyw2au(cmn, epsh, 1) +
          w3stur));
        b4stui = sqmq / sqmh * (sqmq / sqmh - 1.f / 4.f) * (pyw2au(cmn,
          epss, 2) - pyw2au(cmn, epsh, 2) + w3stui);
        b4tusr = sqmq / sqmh * (-2.f / 3.f + (sqmq / sqmh - 1.f /
          4.f) * (pyw2au(cmn, epst, 1) - pyw2au(cmn, epsh, 1) +
          w3tusr));
        b4tusi = sqmq / sqmh * (sqmq / sqmh - 1.f / 4.f) * (pyw2au(cmn,
          epst, 2) - pyw2au(cmn, epsh, 2) + w3tusi);
        b4ustr = sqmq / sqmh * (-2.f / 3.f + (sqmq / sqmh - 1.f /
          4.f) * (pyw2au(cmn, epsu, 1) - pyw2au(cmn, epsh, 1) +
          w3ustr));
        b4usti = sqmq / sqmh * (sqmq / sqmh - 1.f / 4.f) * (pyw2au(cmn,
          epsu, 2) - pyw2au(cmn, epsh, 2) + w3usti);
        a2stur += b2stur + b2sutr;
        a2stui += b2stui + b2suti;
        a2ustr += b2ustr + b2utsr;
        a2usti += b2usti + b2utsi;
        a2tusr += b2tusr + b2tsur;
        a2tusi += b2tusi + b2tsui;
        a4stur += b4stur + b4ustr + b4tusr;
        a4stui += b4stui + b4usti + b4tusi;
        statement_760:;
      }
      facgh = comfac * faca * 3.f / (128.f * fem::pow2(paru(1))) *
        aem / xw * fem::pow3(as) * sqmh / sqmw * fem::pow3(sqmh) / (
        sh * th * uh) * (fem::pow2(a2stur) + fem::pow2(a2stui) +
        fem::pow2(a2ustr) + fem::pow2(a2usti) + fem::pow2(a2tusr) +
        fem::pow2(a2tusi) + fem::pow2(a4stur) + fem::pow2(a4stui));
      facgh = facgh * wids(25, 2);
      if (kfac(1, 21) * kfac(2, 21) == 0) {
        goto statement_770;
      }
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 1;
      sigh(nchn) = facgh;
      statement_770:;
      //C
    }
    else if (isub == 114) {
      //C...g + g -> gamma + gamma.
      asre = 0.f;
      asim = 0.f;
      FEM_DO_SAFE(i, 1, 2 * mstp(1)) {
        ei = kchg(fem::iabs(i), 1) / 3.f;
        sqmq = fem::pow2(pmas(i, 1));
        epss = 4.f * sqmq / sh;
        epst = 4.f * sqmq / th;
        epsu = 4.f * sqmq / uh;
        if (epss + fem::abs(epst) + fem::abs(epsu) < 3.e-6f) {
          a0stur = 1.f + (th - uh) / sh * fem::log(th / uh) + 0.5f * (
            th2 + uh2) / sh2 * (fem::pow2(fem::log(th / uh)) +
            fem::pow2(paru(1)));
          a0stui = 0.f;
          a0tsur = 1.f + (sh - uh) / th * fem::log(-sh / uh) + 0.5f *
            (sh2 + uh2) / th2 * fem::pow2(fem::log(-sh / uh));
          a0tsui = -paru(1) * ((sh - uh) / th + (sh2 + uh2) / th2 *
            fem::log(-sh / uh));
          a0utsr = 1.f + (th - sh) / uh * fem::log(-th / sh) + 0.5f *
            (th2 + sh2) / uh2 * fem::pow2(fem::log(-th / sh));
          a0utsi = paru(1) * ((th - sh) / uh + (th2 + sh2) / uh2 *
            fem::log(-th / sh));
          a1stur = -1.f;
          a1stui = 0.f;
          a2stur = -1.f;
          a2stui = 0.f;
        }
        else {
          bestu = 0.5f * (1.f + fem::sqrt(1.f + epss * th / uh));
          beust = 0.5f * (1.f + fem::sqrt(1.f + epsu * sh / th));
          betus = 0.5f * (1.f + fem::sqrt(1.f + epst * uh / sh));
          beuts = bestu;
          betsu = beust;
          besut = betus;
          a0stur = 1.f + (1.f + 2.f * th / sh) * pyw1au(cmn, epst,
            1) + (1.f + 2.f * uh / sh) * pyw1au(cmn, epsu, 1) +
            0.5f * ((th2 + uh2) / sh2 - epss) * (pyw2au(cmn, epst,
            1) + pyw2au(cmn, epsu, 1)) - 0.25f * epst * (1.f - 0.5f *
            epss) * (pyi3au(cmn, besut, epss, 1) + pyi3au(cmn, besut,
            epst, 1)) - 0.25f * epsu * (1.f - 0.5f * epss) * (pyi3au(cmn,
            bestu, epss, 1) + pyi3au(cmn, bestu, epsu, 1)) + 0.25f * (
            -2.f * (th2 + uh2) / sh2 + 4.f * epss + epst + epsu +
            0.5f * epst * epsu) * (pyi3au(cmn, betsu, epst, 1) + pyi3au(cmn,
            betsu, epsu, 1));
          a0stui = (1.f + 2.f * th / sh) * pyw1au(cmn, epst, 2) + (
            1.f + 2.f * uh / sh) * pyw1au(cmn, epsu, 2) + 0.5f * ((
            th2 + uh2) / sh2 - epss) * (pyw2au(cmn, epst, 2) + pyw2au(cmn,
            epsu, 2)) - 0.25f * epst * (1.f - 0.5f * epss) * (pyi3au(cmn,
            besut, epss, 2) + pyi3au(cmn, besut, epst, 2)) - 0.25f *
            epsu * (1.f - 0.5f * epss) * (pyi3au(cmn, bestu, epss,
            2) + pyi3au(cmn, bestu, epsu, 2)) + 0.25f * (-2.f * (
            th2 + uh2) / sh2 + 4.f * epss + epst + epsu + 0.5f *
            epst * epsu) * (pyi3au(cmn, betsu, epst, 2) + pyi3au(cmn,
            betsu, epsu, 2));
          a0tsur = 1.f + (1.f + 2.f * sh / th) * pyw1au(cmn, epss,
            1) + (1.f + 2.f * uh / th) * pyw1au(cmn, epsu, 1) +
            0.5f * ((sh2 + uh2) / th2 - epst) * (pyw2au(cmn, epss,
            1) + pyw2au(cmn, epsu, 1)) - 0.25f * epss * (1.f - 0.5f *
            epst) * (pyi3au(cmn, betus, epst, 1) + pyi3au(cmn, betus,
            epss, 1)) - 0.25f * epsu * (1.f - 0.5f * epst) * (pyi3au(cmn,
            betsu, epst, 1) + pyi3au(cmn, betsu, epsu, 1)) + 0.25f * (
            -2.f * (sh2 + uh2) / th2 + 4.f * epst + epss + epsu +
            0.5f * epss * epsu) * (pyi3au(cmn, bestu, epss, 1) + pyi3au(cmn,
            bestu, epsu, 1));
          a0tsui = (1.f + 2.f * sh / th) * pyw1au(cmn, epss, 2) + (
            1.f + 2.f * uh / th) * pyw1au(cmn, epsu, 2) + 0.5f * ((
            sh2 + uh2) / th2 - epst) * (pyw2au(cmn, epss, 2) + pyw2au(cmn,
            epsu, 2)) - 0.25f * epss * (1.f - 0.5f * epst) * (pyi3au(cmn,
            betus, epst, 2) + pyi3au(cmn, betus, epss, 2)) - 0.25f *
            epsu * (1.f - 0.5f * epst) * (pyi3au(cmn, betsu, epst,
            2) + pyi3au(cmn, betsu, epsu, 2)) + 0.25f * (-2.f * (
            sh2 + uh2) / th2 + 4.f * epst + epss + epsu + 0.5f *
            epss * epsu) * (pyi3au(cmn, bestu, epss, 2) + pyi3au(cmn,
            bestu, epsu, 2));
          a0utsr = 1.f + (1.f + 2.f * th / uh) * pyw1au(cmn, epst,
            1) + (1.f + 2.f * sh / uh) * pyw1au(cmn, epss, 1) +
            0.5f * ((th2 + sh2) / uh2 - epsu) * (pyw2au(cmn, epst,
            1) + pyw2au(cmn, epss, 1)) - 0.25f * epst * (1.f - 0.5f *
            epsu) * (pyi3au(cmn, beust, epsu, 1) + pyi3au(cmn, beust,
            epst, 1)) - 0.25f * epss * (1.f - 0.5f * epsu) * (pyi3au(cmn,
            beuts, epsu, 1) + pyi3au(cmn, beuts, epss, 1)) + 0.25f * (
            -2.f * (th2 + sh2) / uh2 + 4.f * epsu + epst + epss +
            0.5f * epst * epss) * (pyi3au(cmn, betus, epst, 1) + pyi3au(cmn,
            betus, epss, 1));
          a0utsi = (1.f + 2.f * th / uh) * pyw1au(cmn, epst, 2) + (
            1.f + 2.f * sh / uh) * pyw1au(cmn, epss, 2) + 0.5f * ((
            th2 + sh2) / uh2 - epsu) * (pyw2au(cmn, epst, 2) + pyw2au(cmn,
            epss, 2)) - 0.25f * epst * (1.f - 0.5f * epsu) * (pyi3au(cmn,
            beust, epsu, 2) + pyi3au(cmn, beust, epst, 2)) - 0.25f *
            epss * (1.f - 0.5f * epsu) * (pyi3au(cmn, beuts, epsu,
            2) + pyi3au(cmn, beuts, epss, 2)) + 0.25f * (-2.f * (
            th2 + sh2) / uh2 + 4.f * epsu + epst + epss + 0.5f *
            epst * epss) * (pyi3au(cmn, betus, epst, 2) + pyi3au(cmn,
            betus, epss, 2));
          a1stur = -1.f - 0.25f * (epss + epst + epsu) * (pyw2au(cmn,
            epss, 1) + pyw2au(cmn, epst, 1) + pyw2au(cmn, epsu, 1)) +
            0.25f * (epsu + 0.5f * epss * epst) * (pyi3au(cmn, besut,
            epss, 1) + pyi3au(cmn, besut, epst, 1)) + 0.25f * (epst +
            0.5f * epss * epsu) * (pyi3au(cmn, bestu, epss, 1) + pyi3au(cmn,
            bestu, epsu, 1)) + 0.25f * (epss + 0.5f * epst * epsu) * (
            pyi3au(cmn, betsu, epst, 1) + pyi3au(cmn, betsu, epsu,
            1));
          a1stui = -0.25f * (epss + epst + epsu) * (pyw2au(cmn, epss,
            2) + pyw2au(cmn, epst, 2) + pyw2au(cmn, epsu, 2)) +
            0.25f * (epsu + 0.5f * epss * epst) * (pyi3au(cmn, besut,
            epss, 2) + pyi3au(cmn, besut, epst, 2)) + 0.25f * (epst +
            0.5f * epss * epsu) * (pyi3au(cmn, bestu, epss, 2) + pyi3au(cmn,
            bestu, epsu, 2)) + 0.25f * (epss + 0.5f * epst * epsu) * (
            pyi3au(cmn, betsu, epst, 2) + pyi3au(cmn, betsu, epsu,
            2));
          a2stur = -1.f + 0.125f * epss * epst * (pyi3au(cmn, besut,
            epss, 1) + pyi3au(cmn, besut, epst, 1)) + 0.125f * epss *
            epsu * (pyi3au(cmn, bestu, epss, 1) + pyi3au(cmn, bestu,
            epsu, 1)) + 0.125f * epst * epsu * (pyi3au(cmn, betsu,
            epst, 1) + pyi3au(cmn, betsu, epsu, 1));
          a2stui = 0.125f * epss * epst * (pyi3au(cmn, besut, epss,
            2) + pyi3au(cmn, besut, epst, 2)) + 0.125f * epss *
            epsu * (pyi3au(cmn, bestu, epss, 2) + pyi3au(cmn, bestu,
            epsu, 2)) + 0.125f * epst * epsu * (pyi3au(cmn, betsu,
            epst, 2) + pyi3au(cmn, betsu, epsu, 2));
        }
        asre += fem::pow2(ei) * (a0stur + a0tsur + a0utsr + 4.f *
          a1stur + a2stur);
        asim += fem::pow2(ei) * (a0stui + a0tsui + a0utsi + 4.f *
          a1stui + a2stui);
      }
      facgg = comfac * faca / (8.f * fem::pow2(paru(1))) * fem::pow2(
        as) * fem::pow2(aem) * (fem::pow2(asre) + fem::pow2(asim));
      if (kfac(1, 21) * kfac(2, 21) == 0) {
        goto statement_790;
      }
      nchn++;
      isig(nchn, 1) = 21;
      isig(nchn, 2) = 21;
      isig(nchn, 3) = 1;
      sigh(nchn) = facgg;
      statement_790:;
      //C
    }
    else if (isub == 115) {
      //C...g + g -> gamma + Z0.
      //C
    }
    else if (isub == 116) {
      //C...g + g -> Z0 + Z0.
      //C
    }
    else if (isub == 117) {
      //C...g + g -> W+ + W-.
      //C
    }
    //C
    //C...G: 2 -> 3, tree diagrams.
    //C
  }
  else if (isub <= 140) {
    if (isub == 121) {
      //C...g + g -> f + fb + H0.
      //C
    }
    //C
    //C...H: 2 -> 1, tree diagrams, non-standard model processes.
    //C
  }
  else if (isub <= 160) {
    if (isub == 141) {
      //C...f + fb -> gamma*/Z0/Z'0.
      mint(61) = 2;
      pywidt(cmn, 32, fem::sqrt(sh), wdtp, wdte);
      faczp = comfac * fem::pow2(aem) * 4.f / 9.f;
      FEM_DO_SAFE(i, mina, maxa) {
        if (i == 0 || kfac(1, i) * kfac(2, -i) == 0) {
          goto statement_800;
        }
        ei = kchg(fem::iabs(i), 1) / 3.f;
        ai = fem::sign(1.f, ei);
        vi = ai - 4.f * ei * xw;
        api = fem::sign(1.f, ei);
        vpi = api - 4.f * ei * xw;
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = faczp * (fem::pow2(ei) * vint(111) + ei * vi / (
          8.f * xw * (1.f - xw)) * sh * (sh - sqmz) / (fem::pow2((
          sh - sqmz)) + fem::pow2(gmmz)) * vint(112) + ei * vpi / (
          8.f * xw * (1.f - xw)) * sh * (sh - sqmzp) / (fem::pow2((
          sh - sqmzp)) + fem::pow2(gmmzp)) * vint(113) + (fem::pow2(
          vi) + fem::pow2(ai)) / fem::pow2((16.f * xw * (1.f -
          xw))) * sh2 / (fem::pow2((sh - sqmz)) + fem::pow2(gmmz)) *
          vint(114) + 2.f * (vi * vpi + ai * api) / fem::pow2((16.f *
          xw * (1.f - xw))) * sh2 * ((sh - sqmz) * (sh - sqmzp) +
          gmmz * gmmzp) / ((fem::pow2((sh - sqmz)) + fem::pow2(
          gmmz)) * (fem::pow2((sh - sqmzp)) + fem::pow2(gmmzp))) *
          vint(115) + (fem::pow2(vpi) + fem::pow2(api)) / fem::pow2((
          16.f * xw * (1.f - xw))) * sh2 / (fem::pow2((sh - sqmzp)) +
          fem::pow2(gmmzp)) * vint(116));
        statement_800:;
      }
      //C
    }
    else if (isub == 142) {
      //C...f + fb' -> H+/-.
      pywidt(cmn, 37, fem::sqrt(sh), wdtp, wdte);
      fhc = comfac * fem::pow2((aem / xw)) * 1.f / 48.f * fem::pow2((
        sh / sqmw)) * sh2 / (fem::pow2((sh - sqmhc)) + fem::pow2(
        gmmhc));
      //C'''No construction yet for leptons
      FEM_DO_SAFE(i, 1, mstp(54) / 2) {
        il = 2 * i - 1;
        iu = 2 * i;
        rmql = fem::pow2(pmas(il, 1)) / sh;
        rmqu = fem::pow2(pmas(iu, 1)) / sh;
        fachc = fhc * ((rmql * paru(121) + rmqu / paru(121)) * (1.f -
          rmql - rmqu) - 4.f * rmql * rmqu) / fem::sqrt(fem::max(0.f,
          fem::pow2((1.f - rmql - rmqu)) - 4.f * rmql * rmqu));
        if (kfac(1, il) * kfac(2, -iu) == 0) {
          goto statement_810;
        }
        kchhc = (kchg(il, 1) - kchg(iu, 1)) / 3;
        nchn++;
        isig(nchn, 1) = il;
        isig(nchn, 2) = -iu;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachc * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_810:
        if (kfac(1, -il) * kfac(2, iu) == 0) {
          goto statement_820;
        }
        kchhc = (-kchg(il, 1) + kchg(iu, 1)) / 3;
        nchn++;
        isig(nchn, 1) = -il;
        isig(nchn, 2) = iu;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachc * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_820:
        if (kfac(1, iu) * kfac(2, -il) == 0) {
          goto statement_830;
        }
        kchhc = (kchg(iu, 1) - kchg(il, 1)) / 3;
        nchn++;
        isig(nchn, 1) = iu;
        isig(nchn, 2) = -il;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachc * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_830:
        if (kfac(1, -iu) * kfac(2, il) == 0) {
          goto statement_840;
        }
        kchhc = (-kchg(iu, 1) + kchg(il, 1)) / 3;
        nchn++;
        isig(nchn, 1) = -iu;
        isig(nchn, 2) = il;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachc * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_840:;
      }
      //C
    }
    else if (isub == 143) {
      //C...f + fb -> R.
      pywidt(cmn, 40, fem::sqrt(sh), wdtp, wdte);
      facr = comfac * fem::pow2((aem / xw)) * 1.f / 9.f * sh2 / (
        fem::pow2((sh - sqmr)) + fem::pow2(gmmr));
      FEM_DO_SAFE(i, min1, max1) {
        if (i == 0 || kfac(1, i) == 0) {
          goto statement_860;
        }
        ia = fem::iabs(i);
        FEM_DO_SAFE(j, min2, max2) {
          if (j == 0 || kfac(2, j) == 0) {
            goto statement_850;
          }
          ja = fem::iabs(j);
          if (i * j > 0 || fem::iabs(ia - ja) != 2) {
            goto statement_850;
          }
          nchn++;
          isig(nchn, 1) = i;
          isig(nchn, 2) = j;
          isig(nchn, 3) = 1;
          sigh(nchn) = facr * (wdte(0, 1) + wdte(0, (10 - (i + j)) /
            4) + wdte(0, 4));
          statement_850:;
        }
        statement_860:;
      }
      //C
    }
    //C
    //C...I: 2 -> 2, tree diagrams, non-standard model processes.
    //C
  }
  else {
    if (isub == 161) {
      //C...f + g -> f' + H+/- (q + g -> q' + H+/- only).
      fhcq = comfac * faca * as * aem / xw * 1.f / 24;
      FEM_DO_SAFE(i, 1, mstp(54)) {
        iu = i + fem::mod(i, 2);
        sqmq = fem::pow2(pmas(iu, 1));
        fachcq = fhcq / paru(121) * sqmq / sqmw * (sh / (sqmq - uh) +
          2.f * sqmq * (sqmhc - uh) / fem::pow2((sqmq - uh)) + (
          sqmq - uh) / sh + 2.f * sqmq / (sqmq - uh) + 2.f * (sqmhc -
          uh) / (sqmq - uh) * (sqmhc - sqmq - sh) / sh);
        if (kfac(1, -i) * kfac(2, 21) == 0) {
          goto statement_870;
        }
        kchhc = fem::isign(1, -kchg(i, 1));
        nchn++;
        isig(nchn, 1) = -i;
        isig(nchn, 2) = 21;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachcq * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_870:
        if (kfac(1, i) * kfac(2, 21) == 0) {
          goto statement_880;
        }
        kchhc = fem::isign(1, kchg(i, 1));
        nchn++;
        isig(nchn, 1) = i;
        isig(nchn, 2) = 21;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachcq * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_880:
        if (kfac(1, 21) * kfac(2, -i) == 0) {
          goto statement_890;
        }
        kchhc = fem::isign(1, -kchg(i, 1));
        nchn++;
        isig(nchn, 1) = 21;
        isig(nchn, 2) = -i;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachcq * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_890:
        if (kfac(1, 21) * kfac(2, i) == 0) {
          goto statement_900;
        }
        kchhc = fem::isign(1, kchg(i, 1));
        nchn++;
        isig(nchn, 1) = 21;
        isig(nchn, 2) = i;
        isig(nchn, 3) = 1;
        sigh(nchn) = fachcq * (wdte(0, 1) + wdte(0, (5 - kchhc) / 2) + wdte(0,
          4));
        statement_900:;
      }
      //C
    }
  }
  //C
  //C...Multiply with structure functions.
  if (isub <= 90 || isub >= 96) {
    FEM_DO_SAFE(ichn, 1, nchn) {
      if (mint(41) == 2) {
        kfl1 = isig(ichn, 1);
        if (kfl1 == 21) {
          kfl1 = 0;
        }
        sigh(ichn) = sigh(ichn) * xsfx(1, kfl1);
      }
      if (mint(42) == 2) {
        kfl2 = isig(ichn, 2);
        if (kfl2 == 21) {
          kfl2 = 0;
        }
        sigh(ichn) = sigh(ichn) * xsfx(2, kfl2);
      }
      sigs += sigh(ichn);
    }
  }
  //C
}

struct pymaxi_save
{
  arr<fem::str<4> > cvar;

  pymaxi_save() :
    cvar(dimension(4), fem::fill0)
  {}
};

//C
//C*********************************************************************
//C
void
pymaxi(
  common& cmn)
{
  FEM_CMN_SVE(pymaxi);
  common_write write(cmn);
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int> msub(cmn.msub, dimension(200));
  arr_cref<float> ckin(cmn.ckin, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_cref<int, 2> kfpr(cmn.kfpr, dimension(200, 2));
  arr_ref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  str_arr_cref<1> proc(cmn.proc, dim1(0, 200));
  //
  str_arr_ref<1> cvar(sve.cvar, dimension(4));
  if (is_called_first_time) {
    static const char* values[] = {
      "tau ", "tau'", "y*  ", "cth "
    };
    fem::data_of_type_str(FEM_VALUES_AND_SIZE),
      cvar;
  }
  int isub = fem::int0;
  int istsb = fem::int0;
  int kfr1 = fem::int0;
  float taur1 = fem::float0;
  float gamr1 = fem::float0;
  int kfr2 = fem::int0;
  float taur2 = fem::float0;
  float gamr2 = fem::float0;
  float sqm3 = fem::float0;
  float sqm4 = fem::float0;
  arr_1d<4, int> npts(fem::fill0);
  int ntry = fem::int0;
  int j = fem::int0;
  int mcth = fem::int0;
  int mtaup = fem::int0;
  float cth = fem::float0;
  float taup = fem::float0;
  float sigsam = fem::float0;
  int nacc = fem::int0;
  int itry = fem::int0;
  int mtau = fem::int0;
  int myst = fem::int0;
  arr<int, 2> mvarpt(dimension(200, 4), fem::fill0);
  arr<float, 2> vintpt(dimension(200, 30), fem::fill0);
  int nchn = fem::int0;
  float sigs = fem::float0;
  arr_1d<200, float> sigspt(fem::fill0);
  float taumin = fem::float0;
  float taumax = fem::float0;
  float atau1 = fem::float0;
  float atau2 = fem::float0;
  float atau3 = fem::float0;
  float atau4 = fem::float0;
  float atau5 = fem::float0;
  float atau6 = fem::float0;
  float ystmin = fem::float0;
  float ystmax = fem::float0;
  float ayst0 = fem::float0;
  float ayst1 = fem::float0;
  float ayst3 = fem::float0;
  int ivar = fem::int0;
  int nbin = fem::int0;
  int j1 = fem::int0;
  arr_1d<6, int> narel(fem::fill0);
  arr_1d<6, float> wtrel(fem::fill0);
  arr_1d<6, float> coefu(fem::fill0);
  int j2 = fem::int0;
  arr_2d<6, 6, float> wtmat(fem::fill0);
  int iacc = fem::int0;
  int ibin = fem::int0;
  float tau = fem::float0;
  float taupmn = fem::float0;
  float taupmx = fem::float0;
  float ataup1 = fem::float0;
  float ataup2 = fem::float0;
  float yst = fem::float0;
  float rm34 = fem::float0;
  float rsqm = fem::float0;
  float cthmax = fem::float0;
  float cthmin = fem::float0;
  float acth1 = fem::float0;
  float acth2 = fem::float0;
  float acth3 = fem::float0;
  float acth4 = fem::float0;
  float acth5 = fem::float0;
  int msolv = fem::int0;
  int ired = fem::int0;
  float rqt = fem::float0;
  int icoe = fem::int0;
  float coefsu = fem::float0;
  int ioff = fem::int0;
  arr_1d<4, int> iaccmx(fem::fill0);
  arr_1d<4, float> sigsmx(fem::fill0);
  int nmax = fem::int0;
  int ieq = fem::int0;
  int imv = fem::int0;
  int iin = fem::int0;
  int imax = fem::int0;
  float vtau = fem::float0;
  float vyst = fem::float0;
  float vcth = fem::float0;
  float vtaup = fem::float0;
  int irpt = fem::int0;
  float vvar = fem::float0;
  int mvar = fem::int0;
  float vdel = fem::float0;
  float vmar = fem::float0;
  int imov0 = fem::int0;
  int imov = fem::int0;
  int inew = fem::int0;
  float vnew = fem::float0;
  arr_1d<3, float> sigssm(fem::fill0);
  float sigs11 = fem::float0;
  //C
  //C...Finds optimal set of coefficients for kinematical variable selection
  //C...and the maximum of the part of the differential cross-section used
  //C...in the event weighting.
  //C
  //C...Select subprocess to study: skip cases not applicable.
  vint(143) = 1.f;
  vint(144) = 1.f;
  xsec(0, 1) = 0.f;
  FEM_DO_SAFE(isub, 1, 200) {
    if (isub >= 91 && isub <= 95) {
      xsec(isub, 1) = vint(isub + 11);
      if (msub(isub) != 1) {
        goto statement_350;
      }
      goto statement_340;
    }
    else if (isub == 96) {
      if (mint(43) != 4) {
        goto statement_350;
      }
      if (msub(95) != 1 && mstp(81) <= 0 && mstp(131) <= 0) {
        goto statement_350;
      }
    }
    else if (isub == 11 || isub == 12 || isub == 13 || isub == 28 ||
      isub == 53 || isub == 68) {
      if (msub(isub) != 1 || msub(95) == 1) {
        goto statement_350;
      }
    }
    else {
      if (msub(isub) != 1) {
        goto statement_350;
      }
    }
    mint(1) = isub;
    istsb = iset(isub);
    if (isub == 96) {
      istsb = 2;
    }
    if (mstp(122) >= 2) {
      write(mstu(11),
        "(/,1x,'Coefficient optimization and maximum search for ',"
        "'subprocess no',i4,/,1x,'Coefficient modes     tau',10x,'y*',9x,"
        "'cth',9x,'tau''',7x,'sigma')"),
        isub;
    }
    //C
    //C...Find resonances (explicit or implicit in cross-section).
    mint(72) = 0;
    kfr1 = 0;
    if (istsb == 1 || istsb == 3) {
      kfr1 = kfpr(isub, 1);
    }
    else if (isub >= 71 && isub <= 77) {
      kfr1 = 25;
    }
    if (kfr1 != 0) {
      taur1 = fem::pow2(pmas(kfr1, 1)) / vint(2);
      gamr1 = pmas(kfr1, 1) * pmas(kfr1, 2) / vint(2);
      mint(72) = 1;
      mint(73) = kfr1;
      vint(73) = taur1;
      vint(74) = gamr1;
    }
    if (isub == 141) {
      kfr2 = 23;
      taur2 = fem::pow2(pmas(kfr2, 1)) / vint(2);
      gamr2 = pmas(kfr2, 1) * pmas(kfr2, 2) / vint(2);
      mint(72) = 2;
      mint(74) = kfr2;
      vint(75) = taur2;
      vint(76) = gamr2;
    }
    //C
    //C...Find product masses and minimum pT of process.
    sqm3 = 0.f;
    sqm4 = 0.f;
    mint(71) = 0;
    vint(71) = ckin(3);
    if (istsb == 2 || istsb == 4) {
      if (kfpr(isub, 1) != 0) {
        sqm3 = fem::pow2(pmas(kfpr(isub, 1), 1));
      }
      if (kfpr(isub, 2) != 0) {
        sqm4 = fem::pow2(pmas(kfpr(isub, 2), 1));
      }
      if (fem::min(sqm3, sqm4) < fem::pow2(ckin(6))) {
        mint(71) = 1;
      }
      if (mint(71) == 1) {
        vint(71) = fem::max(ckin(3), ckin(5));
      }
      if (isub == 96 && mstp(82) <= 1) {
        vint(71) = parp(81);
      }
      if (isub == 96 && mstp(82) >= 2) {
        vint(71) = 0.08f * parp(82);
      }
    }
    vint(63) = sqm3;
    vint(64) = sqm4;
    //C
    //C...Number of points for each variable: tau, tau', y*, cos(theta-hat).
    npts(1) = 2 + 2 * mint(72);
    if (mint(43) == 1 && (istsb == 1 || istsb == 2)) {
      npts(1) = 1;
    }
    npts(2) = 1;
    if (mint(43) >= 2 && (istsb == 3 || istsb == 4)) {
      npts(2) = 2;
    }
    npts(3) = 1;
    if (mint(43) == 4) {
      npts(3) = 3;
    }
    npts(4) = 1;
    if (istsb == 2 || istsb == 4) {
      npts(4) = 5;
    }
    ntry = npts(1) * npts(2) * npts(3) * npts(4);
    //C
    //C...Reset coefficients of cross-section weighting.
    FEM_DO_SAFE(j, 1, 20) {
      coef(isub, j) = 0.f;
    }
    coef(isub, 1) = 1.f;
    coef(isub, 7) = 0.5f;
    coef(isub, 8) = 0.5f;
    coef(isub, 10) = 1.f;
    coef(isub, 15) = 1.f;
    mcth = 0;
    mtaup = 0;
    cth = 0.f;
    taup = 0.f;
    sigsam = 0.f;
    //C
    //C...Find limits and select tau, y*, cos(theta-hat) and tau' values,
    //C...in grid of phase space points.
    pyklim(cmn, 1);
    nacc = 0;
    FEM_DO_SAFE(itry, 1, ntry) {
      if (fem::mod(itry - 1, npts(2) * npts(3) * npts(4)) == 0) {
        mtau = 1 + (itry - 1) / (npts(2) * npts(3) * npts(4));
        pykmap(cmn, 1, mtau, 0.5f);
        if (istsb == 3 || istsb == 4) {
          pyklim(cmn, 4);
        }
      }
      if ((istsb == 3 || istsb == 4) && fem::mod(itry - 1, npts(3) *
          npts(4)) == 0) {
        mtaup = 1 + fem::mod((itry - 1) / (npts(3) * npts(4)), npts(2));
        pykmap(cmn, 4, mtaup, 0.5f);
      }
      if (fem::mod(itry - 1, npts(3) * npts(4)) == 0) {
        pyklim(cmn, 2);
      }
      if (fem::mod(itry - 1, npts(4)) == 0) {
        myst = 1 + fem::mod((itry - 1) / npts(4), npts(3));
        pykmap(cmn, 2, myst, 0.5f);
        pyklim(cmn, 3);
      }
      if (istsb == 2 || istsb == 4) {
        mcth = 1 + fem::mod(itry - 1, npts(4));
        pykmap(cmn, 3, mcth, 0.5f);
      }
      if (isub == 96) {
        vint(25) = vint(21) * (1.f - fem::pow2(vint(23)));
      }
      //C
      //C...Calculate and store cross-section.
      mint(51) = 0;
      pyklim(cmn, 0);
      if (mint(51) == 1) {
        goto statement_120;
      }
      nacc++;
      mvarpt(nacc, 1) = mtau;
      mvarpt(nacc, 2) = mtaup;
      mvarpt(nacc, 3) = myst;
      mvarpt(nacc, 4) = mcth;
      FEM_DO_SAFE(j, 1, 30) {
        vintpt(nacc, j) = vint(10 + j);
      }
      pysigh(cmn, nchn, sigs);
      sigspt(nacc) = sigs;
      if (sigs > sigsam) {
        sigsam = sigs;
      }
      if (mstp(122) >= 2) {
        write(mstu(11), "(1x,4i4,f12.8,f12.6,f12.7,f12.8,1p,e12.4)"),
          mtau, mtaup, myst, mcth, vint(21), vint(22), vint(23), vint(26),
          sigs;
      }
      statement_120:;
    }
    if (sigsam == 0.f) {
      write(mstu(11),
        "(1x,'Error: requested subprocess ',i3,' has vanishing ',"
        "'cross-section.',/,1x,'Execution stopped!')"),
        isub;
      FEM_STOP(0);
    }
    //C
    //C...Calculate integrals in tau and y* over maximal phase space limits.
    taumin = vint(11);
    taumax = vint(31);
    atau1 = fem::log(taumax / taumin);
    atau2 = (taumax - taumin) / (taumax * taumin);
    if (npts(1) >= 3) {
      atau3 = fem::log(taumax / taumin * (taumin + taur1) / (taumax +
        taur1)) / taur1;
      atau4 = (fem::atan((taumax - taur1) / gamr1) - fem::atan((
        taumin - taur1) / gamr1)) / gamr1;
    }
    if (npts(1) >= 5) {
      atau5 = fem::log(taumax / taumin * (taumin + taur2) / (taumax +
        taur2)) / taur2;
      atau6 = (fem::atan((taumax - taur2) / gamr2) - fem::atan((
        taumin - taur2) / gamr2)) / gamr2;
    }
    ystmin = 0.5f * fem::log(taumin);
    ystmax = -ystmin;
    ayst0 = ystmax - ystmin;
    ayst1 = 0.5f * fem::pow2((ystmax - ystmin));
    ayst3 = 2.f * (fem::atan(fem::exp(ystmax)) - fem::atan(fem::exp(ystmin)));
    //C
    //C...Reset. Sum up cross-sections in points calculated.
    FEM_DO_SAFE(ivar, 1, 4) {
      if (npts(ivar) == 1) {
        goto statement_230;
      }
      if (isub == 96 && ivar == 4) {
        goto statement_230;
      }
      nbin = npts(ivar);
      FEM_DO_SAFE(j1, 1, nbin) {
        narel(j1) = 0;
        wtrel(j1) = 0.f;
        coefu(j1) = 0.f;
        FEM_DO_SAFE(j2, 1, nbin) {
          wtmat(j1, j2) = 0.f;
        }
      }
      FEM_DO_SAFE(iacc, 1, nacc) {
        ibin = mvarpt(iacc, ivar);
        narel(ibin)++;
        wtrel(ibin) += sigspt(iacc);
        //C
        //C...Sum up tau cross-section pieces in points used.
        if (ivar == 1) {
          tau = vintpt(iacc, 11);
          wtmat(ibin, 1) += 1.f;
          wtmat(ibin, 2) += (atau1 / atau2) / tau;
          if (nbin >= 3) {
            wtmat(ibin, 3) += (atau1 / atau3) / (tau + taur1);
            wtmat(ibin, 4) += (atau1 / atau4) * tau / (fem::pow2((
              tau - taur1)) + fem::pow2(gamr1));
          }
          if (nbin >= 5) {
            wtmat(ibin, 5) += (atau1 / atau5) / (tau + taur2);
            wtmat(ibin, 6) += (atau1 / atau6) * tau / (fem::pow2((
              tau - taur2)) + fem::pow2(gamr2));
          }
          //C
          //C...Sum up tau' cross-section pieces in points used.
        }
        else if (ivar == 2) {
          tau = vintpt(iacc, 11);
          taup = vintpt(iacc, 16);
          taupmn = vintpt(iacc, 6);
          taupmx = vintpt(iacc, 26);
          ataup1 = fem::log(taupmx / taupmn);
          ataup2 = (fem::pow4((1.f - tau / taupmx)) - fem::pow4((
            1.f - tau / taupmn))) / (4.f * tau);
          wtmat(ibin, 1) += 1.f;
          wtmat(ibin, 2) += (ataup1 / ataup2) * fem::pow3((1.f -
            tau / taup)) / taup;
          //C
          //C...Sum up y* and cos(theta-hat) cross-section pieces in points used.
        }
        else if (ivar == 3) {
          yst = vintpt(iacc, 12);
          wtmat(ibin, 1) += (ayst0 / ayst1) * (yst - ystmin);
          wtmat(ibin, 2) += (ayst0 / ayst1) * (ystmax - yst);
          wtmat(ibin, 3) += (ayst0 / ayst3) / fem::cosh(yst);
        }
        else {
          rm34 = 2.f * sqm3 * sqm4 / fem::pow2((vintpt(iacc, 11) * vint(2)));
          rsqm = 1.f + rm34;
          cthmax = fem::sqrt(1.f - 4.f * fem::pow2(vint(71)) / (
            taumax * vint(2)));
          cthmin = -cthmax;
          if (cthmax > 0.9999f) {
            rm34 = fem::max(rm34, 2.f * fem::pow2(vint(71)) / (
              taumax * vint(2)));
          }
          acth1 = cthmax - cthmin;
          acth2 = fem::log(fem::max(rm34, rsqm - cthmin) / fem::max(rm34,
            rsqm - cthmax));
          acth3 = fem::log(fem::max(rm34, rsqm + cthmax) / fem::max(rm34,
            rsqm + cthmin));
          acth4 = 1.f / fem::max(rm34, rsqm - cthmax) - 1.f / fem::max(rm34,
            rsqm - cthmin);
          acth5 = 1.f / fem::max(rm34, rsqm + cthmin) - 1.f / fem::max(rm34,
            rsqm + cthmax);
          cth = vintpt(iacc, 13);
          wtmat(ibin, 1) += 1.f;
          wtmat(ibin, 2) += (acth1 / acth2) / fem::max(rm34, rsqm - cth);
          wtmat(ibin, 3) += (acth1 / acth3) / fem::max(rm34, rsqm + cth);
          wtmat(ibin, 4) += (acth1 / acth4) / fem::pow2(fem::max(rm34,
            rsqm - cth));
          wtmat(ibin, 5) += (acth1 / acth5) / fem::pow2(fem::max(rm34,
            rsqm + cth));
        }
      }
      //C
      //C...Check that equation system solvable; else trivial way out.
      if (mstp(122) >= 2) {
        write(mstu(11),
          "(1x,'Coefficients of equation system to be solved for ',a4)"),
          cvar(ivar);
      }
      msolv = 1;
      FEM_DO_SAFE(ibin, 1, nbin) {
        if (mstp(122) >= 2) {
          {
            write_loop wloop(cmn, mstu(11), "(1x,1p,7e11.3)");
            FEM_DO_SAFE(ired, 1, nbin) {
              wloop, wtmat(ibin, ired);
            }
            wloop, wtrel(ibin);
          }
        }
        if (narel(ibin) == 0) {
          msolv = 0;
        }
      }
      if (msolv == 0) {
        FEM_DO_SAFE(ibin, 1, nbin) {
          coefu(ibin) = 1.f;
        }
        //C
        //C...Solve to find relative importance of cross-section pieces.
      }
      else {
        FEM_DO_SAFE(ired, 1, nbin - 1) {
          FEM_DO_SAFE(ibin, ired + 1, nbin) {
            rqt = wtmat(ibin, ired) / wtmat(ired, ired);
            wtrel(ibin) = wtrel(ibin) - rqt * wtrel(ired);
            FEM_DO_SAFE(icoe, ired, nbin) {
              wtmat(ibin, icoe) = wtmat(ibin, icoe) - rqt * wtmat(ired, icoe);
            }
          }
        }
        FEM_DOSTEP(ired, nbin, 1, -1) {
          FEM_DO_SAFE(icoe, ired + 1, nbin) {
            wtrel(ired) = wtrel(ired) - wtmat(ired, icoe) * coefu(icoe);
          }
          coefu(ired) = wtrel(ired) / wtmat(ired, ired);
        }
      }
      //C
      //C...Normalize coefficients, with piece shared democratically.
      coefsu = 0.f;
      FEM_DO_SAFE(ibin, 1, nbin) {
        coefu(ibin) = fem::max(0.f, coefu(ibin));
        coefsu += coefu(ibin);
      }
      if (ivar == 1) {
        ioff = 0;
      }
      if (ivar == 2) {
        ioff = 14;
      }
      if (ivar == 3) {
        ioff = 6;
      }
      if (ivar == 4) {
        ioff = 9;
      }
      if (coefsu > 0.f) {
        FEM_DO_SAFE(ibin, 1, nbin) {
          coef(isub, ioff + ibin) = parp(121) / nbin + (1.f - parp(
            121)) * coefu(ibin) / coefsu;
        }
      }
      else {
        FEM_DO_SAFE(ibin, 1, nbin) {
          coef(isub, ioff + ibin) = 1.f / nbin;
        }
      }
      if (mstp(122) >= 2) {
        {
          write_loop wloop(cmn, mstu(11), "(1x,'Result for ',a4,':',6f9.4)");
          wloop, cvar(ivar);
          FEM_DO_SAFE(ibin, 1, nbin) {
            wloop, coef(isub, ioff + ibin);
          }
        }
      }
      statement_230:;
    }
    //C
    //C...Find two most promising maxima among points previously determined.
    FEM_DO_SAFE(j, 1, 4) {
      iaccmx(j) = 0;
      sigsmx(j) = 0.f;
    }
    nmax = 0;
    FEM_DO_SAFE(iacc, 1, nacc) {
      FEM_DO_SAFE(j, 1, 30) {
        vint(10 + j) = vintpt(iacc, j);
      }
      pysigh(cmn, nchn, sigs);
      ieq = 0;
      FEM_DO_SAFE(imv, 1, nmax) {
        if (fem::abs(sigs - sigsmx(imv)) < 1e-4f * (sigs + sigsmx(imv))) {
          ieq = imv;
        }
      }
      if (ieq == 0) {
        FEM_DOSTEP(imv, nmax, 1, -1) {
          iin = imv + 1;
          if (sigs <= sigsmx(imv)) {
            goto statement_280;
          }
          iaccmx(imv + 1) = iaccmx(imv);
          sigsmx(imv + 1) = sigsmx(imv);
        }
        iin = 1;
        statement_280:
        iaccmx(iin) = iacc;
        sigsmx(iin) = sigs;
        if (nmax <= 1) {
          nmax++;
        }
      }
    }
    //C
    //C...Read out starting position for search.
    if (mstp(122) >= 2) {
      write(mstu(11),
        "(1x,'Maximum search for given coefficients',/,2x,'MAX VAR ',"
        "'MOD MOV   VNEW',7x,'tau',7x,'y*',8x,'cth',7x,'tau''',7x,'sigma')");
    }
    sigsam = sigsmx(1);
    FEM_DO_SAFE(imax, 1, nmax) {
      iacc = iaccmx(imax);
      mtau = mvarpt(iacc, 1);
      mtaup = mvarpt(iacc, 2);
      myst = mvarpt(iacc, 3);
      mcth = mvarpt(iacc, 4);
      vtau = 0.5f;
      vyst = 0.5f;
      vcth = 0.5f;
      vtaup = 0.5f;
      //C
      //C...Starting point and step size in parameter space.
      FEM_DO_SAFE(irpt, 1, 2) {
        FEM_DO_SAFE(ivar, 1, 4) {
          if (npts(ivar) == 1) {
            goto statement_310;
          }
          if (ivar == 1) {
            vvar = vtau;
          }
          if (ivar == 2) {
            vvar = vtaup;
          }
          if (ivar == 3) {
            vvar = vyst;
          }
          if (ivar == 4) {
            vvar = vcth;
          }
          if (ivar == 1) {
            mvar = mtau;
          }
          if (ivar == 2) {
            mvar = mtaup;
          }
          if (ivar == 3) {
            mvar = myst;
          }
          if (ivar == 4) {
            mvar = mcth;
          }
          if (irpt == 1) {
            vdel = 0.1f;
          }
          if (irpt == 2) {
            vdel = fem::max(0.01f, fem::min(0.05f, vvar - 0.02f, 0.98f - vvar));
          }
          if (irpt == 1) {
            vmar = 0.02f;
          }
          if (irpt == 2) {
            vmar = 0.002f;
          }
          imov0 = 1;
          if (irpt == 1 && ivar == 1) {
            imov0 = 0;
          }
          FEM_DO_SAFE(imov, imov0, 8) {
            //C
            //C...Define new point in parameter space.
            if (imov == 0) {
              inew = 2;
              vnew = vvar;
            }
            else if (imov == 1) {
              inew = 3;
              vnew = vvar + vdel;
            }
            else if (imov == 2) {
              inew = 1;
              vnew = vvar - vdel;
            }
            else if (sigssm(3) >= fem::max(sigssm(1), sigssm(2)) &&
              vvar + 2.f * vdel < 1.f - vmar) {
              vvar += vdel;
              sigssm(1) = sigssm(2);
              sigssm(2) = sigssm(3);
              inew = 3;
              vnew = vvar + vdel;
            }
            else if (sigssm(1) >= fem::max(sigssm(2), sigssm(3)) &&
              vvar - 2.f * vdel > vmar) {
              vvar = vvar - vdel;
              sigssm(3) = sigssm(2);
              sigssm(2) = sigssm(1);
              inew = 1;
              vnew = vvar - vdel;
            }
            else if (sigssm(3) >= sigssm(1)) {
              vdel = 0.5f * vdel;
              vvar += vdel;
              sigssm(1) = sigssm(2);
              inew = 2;
              vnew = vvar;
            }
            else {
              vdel = 0.5f * vdel;
              vvar = vvar - vdel;
              sigssm(3) = sigssm(2);
              inew = 2;
              vnew = vvar;
            }
            //C
            //C...Convert to relevant variables and find derived new limits.
            if (ivar == 1) {
              vtau = vnew;
              pykmap(cmn, 1, mtau, vtau);
              if (istsb == 3 || istsb == 4) {
                pyklim(cmn, 4);
              }
            }
            if (ivar <= 2 && (istsb == 3 || istsb == 4)) {
              if (ivar == 2) {
                vtaup = vnew;
              }
              pykmap(cmn, 4, mtaup, vtaup);
            }
            if (ivar <= 2) {
              pyklim(cmn, 2);
            }
            if (ivar <= 3) {
              if (ivar == 3) {
                vyst = vnew;
              }
              pykmap(cmn, 2, myst, vyst);
              pyklim(cmn, 3);
            }
            if (istsb == 2 || istsb == 4) {
              if (ivar == 4) {
                vcth = vnew;
              }
              pykmap(cmn, 3, mcth, vcth);
            }
            if (isub == 96) {
              vint(25) = vint(21) * (1.f - fem::pow2(vint(23)));
            }
            //C
            //C...Evaluate cross-section. Save new maximum. Final maximum.
            pysigh(cmn, nchn, sigs);
            sigssm(inew) = sigs;
            if (sigs > sigsam) {
              sigsam = sigs;
            }
            if (mstp(122) >= 2) {
              write(mstu(11),
                "(1x,4i4,f8.4,f11.7,f9.3,f11.6,f11.7,1p,e12.4)"),
                imax, ivar, mvar, imov, vnew, vint(21), vint(22),
                vint(23), vint(26), sigs;
            }
          }
          statement_310:;
        }
      }
      if (imax == 1) {
        sigs11 = sigsam;
      }
    }
    xsec(isub, 1) = 1.05f * sigsam;
    statement_340:
    if (isub != 96) {
      xsec(0, 1) += xsec(isub, 1);
    }
    statement_350:;
  }
  //C
  //C...Print summary table.
  if (mstp(122) >= 1) {
    write(mstu(11),
      "(/,1x,8('*'),1x,'PYMAXI: summary of differential ',"
      "'cross-section maximum search',1x,8('*'))");
    write(mstu(11),
      "(/,11x,58('='),/,11x,'I',38x,'I',17x,'I',/,11x,'I  ISUB  ',"
      "'Subprocess name',15x,'I  Maximum value  I',/,11x,'I',38x,'I',17x,'I',"
      "/,11x,58('='),/,11x,'I',38x,'I',17x,'I')");
    FEM_DO_SAFE(isub, 1, 200) {
      if (msub(isub) != 1 && isub != 96) {
        goto statement_360;
      }
      if (isub == 96 && mint(43) != 4) {
        goto statement_360;
      }
      if (isub == 96 && msub(95) != 1 && mstp(81) <= 0) {
        goto statement_360;
      }
      if (msub(95) == 1 && (isub == 11 || isub == 12 || isub == 13 ||
          isub == 28 || isub == 53 || isub == 68)) {
        goto statement_360;
      }
      write(mstu(11),
        "(11x,'I',2x,i3,3x,a28,2x,'I',2x,1p,e12.4,3x,'I')"), isub,
        proc(isub), xsec(isub, 1);
      statement_360:;
    }
    write(mstu(11), "(11x,'I',38x,'I',17x,'I',/,11x,58('='))");
  }
  //C
  //C...Format statements for maximization results.
  //C
}

struct pyovly_save
{
  int imax;
  arr<float> wti;
  float wts;

  pyovly_save() :
    imax(fem::int0),
    wti(dim1(0, 100), fem::fill0),
    wts(fem::float0)
  {}
};

//C
//C*********************************************************************
//C
void
pyovly(
  common& cmn,
  int const& movly)
{
  FEM_CMN_SVE(pyovly);
  common_write write(cmn);
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  //
  int& imax = sve.imax;
  arr_ref<float> wti(sve.wti, dim1(0, 100));
  float& wts = sve.wts;
  float xnave = fem::float0;
  float wtn = fem::float0;
  int i = fem::int0;
  float wtr = fem::float0;
  static const char* format_1000 =
    "(1x,'Warning: requested average number of events per bunch',"
    "'crossing too large, ',1p,e12.4)";
  //C
  //C...Initializes multiplicity distribution and selects mutliplicity
  //C...of overlayed events, i.e. several events occuring at the same
  //C...beam crossing.
  //C
  //C...Sum of allowed cross-sections for overlayed events.
  if (movly == 1) {
    vint(131) = vint(106);
    if (mstp(132) >= 2) {
      vint(131) += vint(104);
    }
    if (mstp(132) >= 3) {
      vint(131) += vint(103);
    }
    if (mstp(132) >= 4) {
      vint(131) += vint(102);
    }
    //C
    //C...Initialize multiplicity distribution for unbiased events.
    if (mstp(133) == 1) {
      xnave = vint(131) * parp(131);
      if (xnave > 40.f) {
        write(mstu(11), format_1000), xnave;
      }
      wti(0) = fem::exp(-fem::min(50.f, xnave));
      wts = 0.f;
      wtn = 0.f;
      FEM_DO_SAFE(i, 1, 100) {
        wti(i) = wti(i - 1) * xnave / i;
        if (i - 2.5f > xnave && wti(i) < 1e-6f) {
          goto statement_110;
        }
        wts += wti(i);
        wtn += wti(i) * i;
        imax = i;
      }
      statement_110:
      vint(132) = xnave;
      vint(133) = wtn / wts;
      vint(134) = wts;
      //C
      //C...Initialize mutiplicity distribution for biased events.
    }
    else if (mstp(133) == 2) {
      xnave = vint(131) * parp(131);
      if (xnave > 40.f) {
        write(mstu(11), format_1000), xnave;
      }
      wti(1) = fem::exp(-fem::min(50.f, xnave)) * xnave;
      wts = wti(1);
      wtn = wti(1);
      FEM_DO_SAFE(i, 2, 100) {
        wti(i) = wti(i - 1) * xnave / (i - 1);
        if (i - 2.5f > xnave && wti(i) < 1e-6f) {
          goto statement_130;
        }
        wts += wti(i);
        wtn += wti(i) * i;
        imax = i;
      }
      statement_130:
      vint(132) = xnave;
      vint(133) = wtn / wts;
      vint(134) = wts;
    }
    //C
    //C...Pick multiplicity of overlayed events.
  }
  else {
    if (mstp(133) == 0) {
      mint(81) = fem::max(1, mstp(134));
    }
    else {
      wtr = wts * rlu(cmn, 0);
      FEM_DO_SAFE(i, 1, imax) {
        mint(81) = i;
        wtr = wtr - wti(i);
        if (wtr <= 0.f) {
          goto statement_150;
        }
      }
      statement_150:;
    }
  }
  //C
  //C...Format statement for error message.
  //C
}

struct pymult_save
{
  int irbin;
  arr<int> nmul;
  float rbin;
  arr<float> sigm;
  float xc2;
  float xt2;
  float xt2fac;
  float xts;

  pymult_save() :
    irbin(fem::int0),
    nmul(dimension(20), fem::fill0),
    rbin(fem::float0),
    sigm(dimension(20), fem::fill0),
    xc2(fem::float0),
    xt2(fem::float0),
    xt2fac(fem::float0),
    xts(fem::float0)
  {}
};

//C
//C*********************************************************************
//C
void
pymult(
  common& cmn,
  int const& mmul)
{
  FEM_CMN_SVE(pymult);
  common_write write(cmn);
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float> ckin(cmn.ckin, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_cref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_ref<int, 2> ngen(cmn.ngen, dim1(0, 200).dim2(3));
  arr_cref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  int& irbin = sve.irbin;
  arr_ref<int> nmul(sve.nmul, dimension(20));
  float& rbin = sve.rbin;
  arr_ref<float> sigm(sve.sigm, dimension(20));
  float& xc2 = sve.xc2;
  float& xt2 = sve.xt2;
  float& xt2fac = sve.xt2fac;
  float& xts = sve.xts;
  int isub = fem::int0;
  float sigsum = fem::float0;
  int ixt2 = fem::int0;
  int itry = fem::int0;
  float rsca = fem::float0;
  float taup = fem::float0;
  float tau = fem::float0;
  float ryst = fem::float0;
  int myst = fem::int0;
  int nchn = fem::int0;
  float sigs = fem::float0;
  float yke = fem::float0;
  float so = fem::float0;
  float xi = fem::float0;
  float yi = fem::float0;
  float xk = fem::float0;
  int iit = fem::int0;
  float xf = fem::float0;
  float yf = fem::float0;
  float sp = fem::float0;
  float sop = fem::float0;
  float deltab = fem::float0;
  float b = fem::float0;
  float ov = fem::float0;
  float cq2 = fem::float0;
  float pacc = fem::float0;
  float yk = fem::float0;
  float rtype = fem::float0;
  float b2 = fem::float0;
  float rncor = fem::float0;
  float sigcor = fem::float0;
  int ibin = fem::int0;
  float sigabv = fem::float0;
  int nmax = fem::int0;
  int nstr = fem::int0;
  int i = fem::int0;
  int kcs = fem::int0;
  int j = fem::int0;
  int ist = fem::int0;
  arr<int, 2> kstr(dimension(500, 2), fem::fill0);
  float x1m = fem::float0;
  float x2m = fem::float0;
  float rflav = fem::float0;
  float pt = fem::float0;
  float phi = fem::float0;
  float cth = fem::float0;
  float dmin = fem::float0;
  int istr = fem::int0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  float dist = fem::float0;
  int ist1 = fem::int0;
  int ist2 = fem::int0;
  int istm = fem::int0;
  //C
  //C...Initializes treatment of multiple interactions, selects kinematics
  //C...of hardest interaction if low-pT physics included in run, and
  //C...generates all non-hardest interactions.
  //C
  //C...Initialization of multiple interaction treatment.
  if (mmul == 1) {
    if (mstp(122) >= 1) {
      write(mstu(11),
        "(/,1x,'****** PYMULT: initialization of multiple inter',"
        "'actions for MSTP(82) =',i2,' ******')"),
        mstp(82);
    }
    isub = 96;
    mint(1) = 96;
    vint(63) = 0.f;
    vint(64) = 0.f;
    vint(143) = 1.f;
    vint(144) = 1.f;
    //C
    //C...Loop over phase space points: xT2 choice in 20 bins.
    statement_100:
    sigsum = 0.f;
    FEM_DO_SAFE(ixt2, 1, 20) {
      nmul(ixt2) = mstp(83);
      sigm(ixt2) = 0.f;
      FEM_DO_SAFE(itry, 1, mstp(83)) {
        rsca = 0.05f * ((21 - ixt2) - rlu(cmn, 0));
        xt2 = vint(149) * (1.f + vint(149)) / (vint(149) + rsca) - vint(149);
        xt2 = fem::max(0.01f * vint(149), xt2);
        vint(25) = xt2;
        //C
        //C...Choose tau and y*. Calculate cos(theta-hat).
        if (rlu(cmn, 0) <= coef(isub, 1)) {
          taup = fem::pow((2.f * (1.f + fem::sqrt(1.f - xt2)) / xt2 - 1.f),
            rlu(cmn, 0));
          tau = xt2 * fem::pow2((1.f + taup)) / (4.f * taup);
        }
        else {
          tau = xt2 * (1.f + fem::pow2(fem::tan(rlu(cmn, 0) *
            fem::atan(fem::sqrt(1.f / xt2 - 1.f)))));
        }
        vint(21) = tau;
        pyklim(cmn, 2);
        ryst = rlu(cmn, 0);
        myst = 1;
        if (ryst > coef(isub, 7)) {
          myst = 2;
        }
        if (ryst > coef(isub, 7) + coef(isub, 8)) {
          myst = 3;
        }
        pykmap(cmn, 2, myst, rlu(cmn, 0));
        vint(23) = fem::sqrt(fem::max(0.f, 1.f - xt2 / tau)) * fem::pow((-1),
          fem::fint(1.5f + rlu(cmn, 0)));
        //C
        //C...Calculate differential cross-section.
        vint(71) = 0.5f * vint(1) * fem::sqrt(xt2);
        pysigh(cmn, nchn, sigs);
        sigm(ixt2) += sigs;
      }
      sigsum += sigm(ixt2);
    }
    sigsum = sigsum / (20.f * mstp(83));
    //C
    //C...Reject result if sigma(parton-parton) is smaller than hadronic one.
    if (sigsum < 1.1f * vint(106)) {
      if (mstp(122) >= 1) {
        write(mstu(11),
          "(8x,'pT0 =',f5.2,' GeV gives sigma(parton-parton) =',1p,e9.2,"
          "' mb: rejected')"),
          parp(82), sigsum;
      }
      parp(82) = 0.9f * parp(82);
      vint(149) = 4.f * fem::pow2(parp(82)) / vint(2);
      goto statement_100;
    }
    if (mstp(122) >= 1) {
      write(mstu(11),
        "(8x,'pT0 =',f5.2,' GeV gives sigma(parton-parton) =',1p,e9.2,"
        "' mb: accepted')"),
        parp(82), sigsum;
    }
    //C
    //C...Start iteration to find k factor.
    yke = sigsum / vint(106);
    so = 0.5f;
    xi = 0.f;
    yi = 0.f;
    xk = 0.5f;
    iit = 0;
    statement_130:
    if (iit == 0) {
      xk = 2.f * xk;
    }
    else if (iit == 1) {
      xk = 0.5f * xk;
    }
    else {
      xk = xi + (yke - yi) * (xf - xi) / (yf - yi);
    }
    //C
    //C...Evaluate overlap integrals.
    if (mstp(82) == 2) {
      sp = 0.5f * paru(1) * (1.f - fem::exp(-xk));
      sop = sp / paru(1);
    }
    else {
      if (mstp(82) == 3) {
        deltab = 0.02f;
      }
      if (mstp(82) == 4) {
        deltab = fem::min(0.01f, 0.05f * parp(84));
      }
      sp = 0.f;
      sop = 0.f;
      b = -0.5f * deltab;
      statement_140:
      b += deltab;
      if (mstp(82) == 3) {
        ov = fem::exp(-fem::pow2(b)) / paru(2);
      }
      else {
        cq2 = fem::pow2(parp(84));
        ov = (fem::pow2((1.f - parp(83))) * fem::exp(-fem::min(100.f,
          fem::pow2(b))) + 2.f * parp(83) * (1.f - parp(83)) * 2.f / (
          1.f + cq2) * fem::exp(-fem::min(100.f, fem::pow2(b) * 2.f /
          (1.f + cq2))) + fem::pow2(parp(83)) / cq2 * fem::exp(-fem::min(100.f,
          fem::pow2(b) / cq2))) / paru(2);
      }
      pacc = 1.f - fem::exp(-fem::min(100.f, paru(1) * xk * ov));
      sp += paru(2) * b * deltab * pacc;
      sop += paru(2) * b * deltab * ov * pacc;
      if (b < 1.f || b * pacc > 1e-6f) {
        goto statement_140;
      }
    }
    yk = paru(1) * xk * so / sp;
    //C
    //C...Continue iteration until convergence.
    if (yk < yke) {
      xi = xk;
      yi = yk;
      if (iit == 1) {
        iit = 2;
      }
    }
    else {
      xf = xk;
      yf = yk;
      if (iit == 0) {
        iit = 1;
      }
    }
    if (fem::abs(yk - yke) >= 1e-5f * yke) {
      goto statement_130;
    }
    //C
    //C...Store some results for subsequent use.
    vint(145) = sigsum;
    vint(146) = sop / so;
    vint(147) = sop / sp;
    //C
    //C...Initialize iteration in xT2 for hardest interaction.
  }
  else if (mmul == 2) {
    if (mstp(82) <= 0) {
    }
    else if (mstp(82) == 1) {
      xt2 = 1.f;
      xt2fac = xsec(96, 1) / vint(106) * vint(149) / (1.f - vint(149));
    }
    else if (mstp(82) == 2) {
      xt2 = 1.f;
      xt2fac = vint(146) * xsec(96, 1) / vint(106) * vint(149) * (
        1.f + vint(149));
    }
    else {
      xc2 = 4.f * fem::pow2(ckin(3)) / vint(2);
      if (ckin(3) <= ckin(5) || mint(82) >= 2) {
        xc2 = 0.f;
      }
    }
    //C
  }
  else if (mmul == 3) {
    //C...Low-pT or multiple interactions (first semihard interaction):
    //C...choose xT2 according to dpT2/pT2**2*exp(-(sigma above pT2)/norm)
    //C...or (MSTP(82)>=2) dpT2/(pT2+pT0**2)**2*exp(-....).
    isub = mint(1);
    if (mstp(82) <= 0) {
      xt2 = 0.f;
    }
    else if (mstp(82) == 1) {
      xt2 = xt2fac * xt2 / (xt2fac - xt2 * fem::log(rlu(cmn, 0)));
    }
    else if (mstp(82) == 2) {
      if (xt2 < 1.f && fem::exp(-xt2fac * xt2 / (vint(149) * (xt2 +
          vint(149)))) > rlu(cmn, 0)) {
        xt2 = 1.f;
      }
      if (xt2 >= 1.f) {
        xt2 = (1.f + vint(149)) * xt2fac / (xt2fac - (1.f + vint(
          149)) * fem::log(1.f - rlu(cmn, 0) * (1.f - fem::exp(
          -xt2fac / (vint(149) * (1.f + vint(149))))))) - vint(149);
      }
      else {
        xt2 = -xt2fac / fem::log(fem::exp(-xt2fac / (xt2 + vint(
          149))) + rlu(cmn, 0) * (fem::exp(-xt2fac / vint(149)) -
          fem::exp(-xt2fac / (xt2 + vint(149))))) - vint(149);
      }
      xt2 = fem::max(0.01f * vint(149), xt2);
    }
    else {
      xt2 = (xc2 + vint(149)) * (1.f + vint(149)) / (1.f + vint(149) - rlu(cmn,
        0) * (1.f - xc2)) - vint(149);
      xt2 = fem::max(0.01f * vint(149), xt2);
    }
    vint(25) = xt2;
    //C
    //C...Low-pT: choose xT2, tau, y* and cos(theta-hat) fixed.
    if (mstp(82) <= 1 && xt2 < vint(149)) {
      if (mint(82) == 1) {
        ngen(0, 1) = ngen(0, 1) - 1;
      }
      if (mint(82) == 1) {
        ngen(isub, 1) = ngen(isub, 1) - 1;
      }
      isub = 95;
      mint(1) = isub;
      vint(21) = 0.01f * vint(149);
      vint(22) = 0.f;
      vint(23) = 0.f;
      vint(25) = 0.01f * vint(149);
      //C
    }
    else {
      //C...Multiple interactions (first semihard interaction).
      //C...Choose tau and y*. Calculate cos(theta-hat).
      if (rlu(cmn, 0) <= coef(isub, 1)) {
        taup = fem::pow((2.f * (1.f + fem::sqrt(1.f - xt2)) / xt2 - 1.f),
          rlu(cmn, 0));
        tau = xt2 * fem::pow2((1.f + taup)) / (4.f * taup);
      }
      else {
        tau = xt2 * (1.f + fem::pow2(fem::tan(rlu(cmn, 0) * fem::atan(
          fem::sqrt(1.f / xt2 - 1.f)))));
      }
      vint(21) = tau;
      pyklim(cmn, 2);
      ryst = rlu(cmn, 0);
      myst = 1;
      if (ryst > coef(isub, 7)) {
        myst = 2;
      }
      if (ryst > coef(isub, 7) + coef(isub, 8)) {
        myst = 3;
      }
      pykmap(cmn, 2, myst, rlu(cmn, 0));
      vint(23) = fem::sqrt(fem::max(0.f, 1.f - xt2 / tau)) * fem::pow((-1),
        fem::fint(1.5f + rlu(cmn, 0)));
    }
    vint(71) = 0.5f * vint(1) * fem::sqrt(vint(25));
    //C
    //C...Store results of cross-section calculation.
  }
  else if (mmul == 4) {
    isub = mint(1);
    xts = vint(25);
    if (iset(isub) == 1) {
      xts = vint(21);
    }
    if (iset(isub) == 2) {
      xts = (4.f * vint(48) + 2.f * vint(63) + 2.f * vint(64)) / vint(2);
    }
    if (iset(isub) == 3 || iset(isub) == 4) {
      xts = vint(26);
    }
    rbin = fem::max(0.000001f, fem::min(0.999999f, xts * (1.f + vint(
      149)) / (xts + vint(149))));
    irbin = fem::fint(1.f + 20.f * rbin);
    if (isub == 96) {
      nmul(irbin)++;
    }
    if (isub == 96) {
      sigm(irbin) += vint(153);
    }
    //C
    //C...Choose impact parameter.
  }
  else if (mmul == 5) {
    if (mstp(82) == 3) {
      vint(148) = rlu(cmn, 0) / (paru(2) * vint(147));
    }
    else {
      rtype = rlu(cmn, 0);
      cq2 = fem::pow2(parp(84));
      if (rtype < fem::pow2((1.f - parp(83)))) {
        b2 = -fem::log(rlu(cmn, 0));
      }
      else if (rtype < 1.f - fem::pow2(parp(83))) {
        b2 = -0.5f * (1.f + cq2) * fem::log(rlu(cmn, 0));
      }
      else {
        b2 = -cq2 * fem::log(rlu(cmn, 0));
      }
      vint(148) = (fem::pow2((1.f - parp(83))) * fem::exp(-fem::min(100.f,
        b2)) + 2.f * parp(83) * (1.f - parp(83)) * 2.f / (1.f +
        cq2) * fem::exp(-fem::min(100.f, b2 * 2.f / (1.f + cq2))) +
        fem::pow2(parp(83)) / cq2 * fem::exp(-fem::min(100.f, b2 /
        cq2))) / (paru(2) * vint(147));
    }
    //C
    //C...Multiple interactions (variable impact parameter) : reject with
    //C...probability exp(-overlap*cross-section above pT/normalization).
    rncor = (irbin - 20.f * rbin) * nmul(irbin);
    sigcor = (irbin - 20.f * rbin) * sigm(irbin);
    FEM_DO_SAFE(ibin, irbin + 1, 20) {
      rncor += nmul(ibin);
      sigcor += sigm(ibin);
    }
    sigabv = (sigcor / rncor) * vint(149) * (1.f - xts) / (xts + vint(149));
    vint(150) = fem::exp(-fem::min(100.f, vint(146) * vint(148) *
      sigabv / vint(106)));
    //C
    //C...Generate additional multiple semihard interactions.
  }
  else if (mmul == 6) {
    //C
    //C...Reconstruct strings in hard scattering.
    isub = mint(1);
    nmax = mint(84) + 4;
    if (iset(isub) == 1) {
      nmax = mint(84) + 2;
    }
    nstr = 0;
    FEM_DO_SAFE(i, mint(84) + 1, nmax) {
      kcs = kchg(lucomp(cmn, k(i, 2)), 2) * fem::isign(1, k(i, 2));
      if (kcs == 0) {
        goto statement_170;
      }
      FEM_DO_SAFE(j, 1, 4) {
        if (kcs == 1 && (j == 2 || j == 4)) {
          goto statement_160;
        }
        if (kcs ==  - 1 && (j == 1 || j == 3)) {
          goto statement_160;
        }
        if (j <= 2) {
          ist = fem::mod(k(i, j + 3) / mstu(5), mstu(5));
        }
        else {
          ist = fem::mod(k(i, j + 1), mstu(5));
        }
        if (ist < mint(84) || ist > i) {
          goto statement_160;
        }
        if (kchg(lucomp(cmn, k(ist, 2)), 2) == 0) {
          goto statement_160;
        }
        nstr++;
        if (j == 1 || j == 4) {
          kstr(nstr, 1) = i;
          kstr(nstr, 2) = ist;
        }
        else {
          kstr(nstr, 1) = ist;
          kstr(nstr, 2) = i;
        }
        statement_160:;
      }
      statement_170:;
    }
    //C
    //C...Set up starting values for iteration in xT2.
    xt2 = vint(25);
    if (iset(isub) == 1) {
      xt2 = vint(21);
    }
    if (iset(isub) == 2) {
      xt2 = (4.f * vint(48) + 2.f * vint(63) + 2.f * vint(64)) / vint(2);
    }
    if (iset(isub) == 3 || iset(isub) == 4) {
      xt2 = vint(26);
    }
    isub = 96;
    mint(1) = 96;
    if (mstp(82) <= 1) {
      xt2fac = xsec(isub, 1) * vint(149) / ((1.f - vint(149)) * vint(106));
    }
    else {
      xt2fac = vint(146) * vint(148) * xsec(isub, 1) / vint(106) *
        vint(149) * (1.f + vint(149));
    }
    vint(63) = 0.f;
    vint(64) = 0.f;
    vint(151) = 0.f;
    vint(152) = 0.f;
    vint(143) = 1.f - vint(141);
    vint(144) = 1.f - vint(142);
    //C
    //C...Iterate downwards in xT2.
    statement_180:
    if (mstp(82) <= 1) {
      xt2 = xt2fac * xt2 / (xt2fac - xt2 * fem::log(rlu(cmn, 0)));
      if (xt2 < vint(149)) {
        goto statement_220;
      }
    }
    else {
      if (xt2 <= 0.01f * vint(149)) {
        goto statement_220;
      }
      xt2 = xt2fac * (xt2 + vint(149)) / (xt2fac - (xt2 + vint(
        149)) * fem::log(rlu(cmn, 0))) - vint(149);
      if (xt2 <= 0.f) {
        goto statement_220;
      }
      xt2 = fem::max(0.01f * vint(149), xt2);
    }
    vint(25) = xt2;
    //C
    //C...Choose tau and y*. Calculate cos(theta-hat).
    if (rlu(cmn, 0) <= coef(isub, 1)) {
      taup = fem::pow((2.f * (1.f + fem::sqrt(1.f - xt2)) / xt2 - 1.f),
        rlu(cmn, 0));
      tau = xt2 * fem::pow2((1.f + taup)) / (4.f * taup);
    }
    else {
      tau = xt2 * (1.f + fem::pow2(fem::tan(rlu(cmn, 0) * fem::atan(
        fem::sqrt(1.f / xt2 - 1.f)))));
    }
    vint(21) = tau;
    pyklim(cmn, 2);
    ryst = rlu(cmn, 0);
    myst = 1;
    if (ryst > coef(isub, 7)) {
      myst = 2;
    }
    if (ryst > coef(isub, 7) + coef(isub, 8)) {
      myst = 3;
    }
    pykmap(cmn, 2, myst, rlu(cmn, 0));
    vint(23) = fem::sqrt(fem::max(0.f, 1.f - xt2 / tau)) * fem::pow((-1),
      fem::fint(1.5f + rlu(cmn, 0)));
    //C
    //C...Check that x not used up. Accept or reject kinematical variables.
    x1m = fem::sqrt(tau) * fem::exp(vint(22));
    x2m = fem::sqrt(tau) * fem::exp(-vint(22));
    if (vint(143) - x1m < 0.01f || vint(144) - x2m < 0.01f) {
      goto statement_180;
    }
    vint(71) = 0.5f * vint(1) * fem::sqrt(xt2);
    pysigh(cmn, nchn, sigs);
    if (sigs < xsec(isub, 1) * rlu(cmn, 0)) {
      goto statement_180;
    }
    //C
    //C...Reset K, P and V vectors. Select some variables.
    FEM_DO_SAFE(i, n + 1, n + 2) {
      FEM_DO_SAFE(j, 1, 5) {
        k(i, j) = 0;
        p(i, j) = 0.f;
        v(i, j) = 0.f;
      }
    }
    rflav = rlu(cmn, 0);
    pt = 0.5f * vint(1) * fem::sqrt(xt2);
    phi = paru(2) * rlu(cmn, 0);
    cth = vint(23);
    //C
    //C...Add first parton to event record.
    k(n + 1, 1) = 3;
    k(n + 1, 2) = 21;
    if (rflav >= fem::max(parp(85), parp(86))) {
      k(n + 1, 2) = 1 + fem::fint((2.f + parj(2)) * rlu(cmn, 0));
    }
    p(n + 1, 1) = pt * fem::cos(phi);
    p(n + 1, 2) = pt * fem::sin(phi);
    p(n + 1, 3) = 0.25f * vint(1) * (vint(41) * (1.f + cth) - vint(
      42) * (1.f - cth));
    p(n + 1, 4) = 0.25f * vint(1) * (vint(41) * (1.f + cth) + vint(
      42) * (1.f - cth));
    p(n + 1, 5) = 0.f;
    //C
    //C...Add second parton to event record.
    k(n + 2, 1) = 3;
    k(n + 2, 2) = 21;
    if (k(n + 1, 2) != 21) {
      k(n + 2, 2) = -k(n + 1, 2);
    }
    p(n + 2, 1) = -p(n + 1, 1);
    p(n + 2, 2) = -p(n + 1, 2);
    p(n + 2, 3) = 0.25f * vint(1) * (vint(41) * (1.f - cth) - vint(
      42) * (1.f + cth));
    p(n + 2, 4) = 0.25f * vint(1) * (vint(41) * (1.f - cth) + vint(
      42) * (1.f + cth));
    p(n + 2, 5) = 0.f;
    //C
    if (rflav < parp(85) && nstr >= 1) {
      //C....Choose relevant string pieces to place gluons on.
      FEM_DO_SAFE(i, n + 1, n + 2) {
        dmin = 1e8f;
        FEM_DO_SAFE(istr, 1, nstr) {
          i1 = kstr(istr, 1);
          i2 = kstr(istr, 2);
          dist = (p(i, 4) * p(i1, 4) - p(i, 1) * p(i1, 1) - p(i, 2) * p(i1,
            2) - p(i, 3) * p(i1, 3)) * (p(i, 4) * p(i2, 4) - p(i,
            1) * p(i2, 1) - p(i, 2) * p(i2, 2) - p(i, 3) * p(i2,
            3)) / fem::max(1.f, p(i1, 4) * p(i2, 4) - p(i1, 1) * p(i2,
            1) - p(i1, 2) * p(i2, 2) - p(i1, 3) * p(i2, 3));
          if (istr == 1 || dist < dmin) {
            dmin = dist;
            ist1 = i1;
            ist2 = i2;
            istm = istr;
          }
        }
        //C
        //C....Colour flow adjustments, new string pieces.
        if (k(ist1, 4) / mstu(5) == ist2) {
          k(ist1, 4) = mstu(5) * i + fem::mod(k(ist1, 4), mstu(5));
        }
        if (fem::mod(k(ist1, 5), mstu(5)) == ist2) {
          k(ist1, 5) = mstu(5) * (k(ist1, 5) / mstu(5)) + i;
        }
        k(i, 5) = mstu(5) * ist1;
        k(i, 4) = mstu(5) * ist2;
        if (k(ist2, 5) / mstu(5) == ist1) {
          k(ist2, 5) = mstu(5) * i + fem::mod(k(ist2, 5), mstu(5));
        }
        if (fem::mod(k(ist2, 4), mstu(5)) == ist1) {
          k(ist2, 4) = mstu(5) * (k(ist2, 4) / mstu(5)) + i;
        }
        kstr(istm, 2) = i;
        kstr(nstr + 1, 1) = i;
        kstr(nstr + 1, 2) = ist2;
        nstr++;
      }
      //C
      //C...String drawing and colour flow for gluon loop.
    }
    else if (k(n + 1, 2) == 21) {
      k(n + 1, 4) = mstu(5) * (n + 2);
      k(n + 1, 5) = mstu(5) * (n + 2);
      k(n + 2, 4) = mstu(5) * (n + 1);
      k(n + 2, 5) = mstu(5) * (n + 1);
      kstr(nstr + 1, 1) = n + 1;
      kstr(nstr + 1, 2) = n + 2;
      kstr(nstr + 2, 1) = n + 2;
      kstr(nstr + 2, 2) = n + 1;
      nstr += 2;
      //C
      //C...String drawing and colour flow for q-qbar pair.
    }
    else {
      k(n + 1, 4) = mstu(5) * (n + 2);
      k(n + 2, 5) = mstu(5) * (n + 1);
      kstr(nstr + 1, 1) = n + 1;
      kstr(nstr + 1, 2) = n + 2;
      nstr++;
    }
    //C
    //C...Update remaining energy; iterate.
    n += 2;
    if (n > mstu(4) - mstu(32) - 10) {
      luerrm(cmn, 11, "(PYMULT:) no more memory left in LUJETS");
      if (mstu(21) >= 1) {
        return;
      }
    }
    mint(31)++;
    vint(151) += vint(41);
    vint(152) += vint(42);
    vint(143) = vint(143) - vint(41);
    vint(144) = vint(144) - vint(42);
    if (mint(31) < 240) {
      goto statement_180;
    }
    statement_220:;
  }
  //C
  //C...Format statements for printout.
  //C
}

struct pyinit_save
{
  arr<fem::str<6> > chlh;
  arr<fem::str<3> > chmo;

  pyinit_save() :
    chlh(dimension(2), fem::fill0),
    chmo(dimension(12), fem::fill0)
  {}
};

void
pyinit(
  common& cmn,
  str_cref frame,
  str_cref beam,
  str_cref target,
  float const& win)
{
  FEM_CMN_SVE(pyinit);
  common_write write(cmn);
  // COMMON ludat1
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<float> paru(cmn.paru, dimension(200));
  arr_ref<float> parj(cmn.parj, dimension(200));
  // COMMON ludat2
  arr_cref<float, 2> vckm(cmn.vckm, dimension(4, 4));
  // COMMON ludat3
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_ref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  // COMMON pysubs
  int& msel = cmn.msel;
  arr_ref<int> msub(cmn.msub, dimension(200));
  arr_cref<float> ckin(cmn.ckin, dimension(200));
  // COMMON pypars
  arr_ref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<float> parp(cmn.parp, dimension(200));
  // COMMON pyint1
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  // COMMON pyint2
  arr_cref<int> iset(cmn.iset, dimension(200));
  // COMMON pyint5
  arr_ref<int, 2> ngen(cmn.ngen, dim1(0, 200).dim2(3));
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  // SAVE
  str_arr_ref<1> chlh(sve.chlh, dimension(2));
  //
  str_arr_ref<1> chmo(sve.chmo, dimension(12));
  if (is_called_first_time) {
    {
      static const char* values[] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
          "Sep", "Oct", "Nov", "Dec"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chmo;
    }
    {
      static const char* values[] = {
        "lepton", "hadron"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        chlh;
    }
  }
  //C
  //C...Initializes the generation procedure; finds maxima of the
  //C...differential cross-sections to be used for weighting.
  //C
  //Clin-12/2012 correct NN differential cross section in HIJING:
  write(mstu(11), star), "In PYINIT: BEAM,TARGET= ", beam, target;
  //C
  //C...Write headers.
  //C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1000) MSTP(181),MSTP(182),
  //C     &MSTP(185),CHMO(MSTP(184)),MSTP(183)
  lulist(cmn, 0);
  //C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1100)
  //C
  //C...Identify beam and target particles and initialize kinematics.
  fem::str<8> chfram = frame + str_cref(" ");
  fem::str<8> chbeam = beam + str_cref(" ");
  fem::str<8> chtarg = target + str_cref(" ");
  pyinki(cmn, chfram, chbeam, chtarg, win);
  //C
  //C...Select partonic subprocesses to be included in the simulation.
  int i = fem::int0;
  if (msel != 0) {
    FEM_DO_SAFE(i, 1, 200) {
      msub(i) = 0;
    }
  }
  int j = fem::int0;
  if (mint(43) == 1 && (msel == 1 || msel == 2)) {
    //C...Lepton+lepton -> gamma/Z0 or W.
    if (mint(11) + mint(12) == 0) {
      msub(1) = 1;
    }
    if (mint(11) + mint(12) != 0) {
      msub(2) = 1;
    }
  }
  else if (msel == 1) {
    //C...High-pT QCD processes:
    msub(11) = 1;
    msub(12) = 1;
    msub(13) = 1;
    msub(28) = 1;
    msub(53) = 1;
    msub(68) = 1;
    if (mstp(82) <= 1 && ckin(3) < parp(81)) {
      msub(95) = 1;
    }
    if (mstp(82) >= 2 && ckin(3) < parp(82)) {
      msub(95) = 1;
    }
  }
  else if (msel == 2) {
    //C...All QCD processes:
    msub(11) = 1;
    msub(12) = 1;
    msub(13) = 1;
    msub(28) = 1;
    msub(53) = 1;
    msub(68) = 1;
    msub(91) = 1;
    msub(92) = 1;
    msub(93) = 1;
    msub(95) = 1;
  }
  else if (msel >= 4 && msel <= 8) {
    //C...Heavy quark production.
    msub(81) = 1;
    msub(82) = 1;
    FEM_DO_SAFE(j, 1, fem::min(8, mdcy(21, 3))) {
      mdme(mdcy(21, 2) + j - 1, 1) = 0;
    }
    mdme(mdcy(21, 2) + msel - 1, 1) = 1;
  }
  else if (msel == 10) {
    //C...Prompt photon production:
    msub(14) = 1;
    msub(18) = 1;
    msub(29) = 1;
  }
  else if (msel == 11) {
    //C...Z0/gamma* production:
    msub(1) = 1;
  }
  else if (msel == 12) {
    //C...W+/- production:
    msub(2) = 1;
  }
  else if (msel == 13) {
    //C...Z0 + jet:
    msub(15) = 1;
    msub(30) = 1;
  }
  else if (msel == 14) {
    //C...W+/- + jet:
    msub(16) = 1;
    msub(31) = 1;
  }
  else if (msel == 15) {
    //C...Z0 & W+/- pair production:
    msub(19) = 1;
    msub(20) = 1;
    msub(22) = 1;
    msub(23) = 1;
    msub(25) = 1;
  }
  else if (msel == 16) {
    //C...H0 production:
    msub(3) = 1;
    msub(5) = 1;
    msub(8) = 1;
    msub(102) = 1;
  }
  else if (msel == 17) {
    //C...H0 & Z0 or W+/- pair production:
    msub(24) = 1;
    msub(26) = 1;
  }
  else if (msel == 21) {
    //C...Z'0 production:
    msub(141) = 1;
  }
  else if (msel == 22) {
    //C...H+/- production:
    msub(142) = 1;
  }
  else if (msel == 23) {
    //C...R production:
    msub(143) = 1;
  }
  //C
  //C...Count number of subprocesses on.
  mint(44) = 0;
  int isub = fem::int0;
  FEM_DO_SAFE(isub, 1, 200) {
    if (mint(43) < 4 && isub >= 91 && isub <= 96 && msub(isub) == 1) {
      write(mstu(11),
        "(1x,'Error: process number ',i3,' not meaningful for ',a6,'-',a6,"
        "' interactions.',/,1x,'Execution stopped!')"),
        isub, chlh(mint(41)), chlh(mint(42));
      FEM_STOP(0);
    }
    else if (msub(isub) == 1 && iset(isub) ==  - 1) {
      write(mstu(11),
        "(1x,'Error: requested subprocess',i4,' not implemented.',/,1x,"
        "'Execution stopped!')"),
        isub;
      FEM_STOP(0);
    }
    else if (msub(isub) == 1 && iset(isub) <=  - 2) {
      write(mstu(11),
        "(1x,'Error: requested subprocess',i4,' not existing.',/,1x,"
        "'Execution stopped!')"),
        isub;
      FEM_STOP(0);
    }
    else if (msub(isub) == 1) {
      mint(44)++;
    }
  }
  if (mint(44) == 0) {
    write(mstu(11),
      "(1x,'Error: no subprocess switched on.',/,1x,'Execution stopped.')");
    FEM_STOP(0);
  }
  mint(45) = mint(44) - msub(91) - msub(92) - msub(93) - msub(94);
  //C
  //C...Maximum 4 generations; set maximum number of allowed flavours.
  mstp(1) = fem::min(4, mstp(1));
  mstu(114) = fem::min(mstu(114), 2 * mstp(1));
  mstp(54) = fem::min(mstp(54), 2 * mstp(1));
  //C
  //C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton.
  int ia = fem::int0;
  int ib = fem::int0;
  int ipm = fem::int0;
  int idc = fem::int0;
  FEM_DO_SAFE(i, -20, 20) {
    vint(180 + i) = 0.f;
    ia = fem::iabs(i);
    if (ia >= 1 && ia <= 2 * mstp(1)) {
      FEM_DO_SAFE(j, 1, mstp(1)) {
        ib = 2 * j - 1 + fem::mod(ia, 2);
        ipm = (5 - fem::isign(1, i)) / 2;
        idc = j + mdcy(ia, 2) + 2;
        if (mdme(idc, 1) == 1 || mdme(idc, 1) == ipm) {
          vint(180 + i) += vckm((ia + 1) / 2, (ib + 1) / 2);
        }
      }
    }
    else if (ia >= 11 && ia <= 10 + 2 * mstp(1)) {
      vint(180 + i) = 1.f;
    }
  }
  //C
  //C...Choose Lambda value to use in alpha-strong.
  mstu(111) = mstp(2);
  float alam = fem::float0;
  if (mstp(3) >= 1) {
    alam = parp(1);
    if (mstp(51) == 1) {
      alam = 0.2f;
    }
    if (mstp(51) == 2) {
      alam = 0.29f;
    }
    if (mstp(51) == 3) {
      alam = 0.2f;
    }
    if (mstp(51) == 4) {
      alam = 0.4f;
    }
    if (mstp(51) == 11) {
      alam = 0.16f;
    }
    if (mstp(51) == 12) {
      alam = 0.26f;
    }
    if (mstp(51) == 13) {
      alam = 0.36f;
    }
    parp(1) = alam;
    parp(61) = alam;
    paru(112) = alam;
    parj(81) = alam;
  }
  //C
  //C...Initialize widths and partial widths for resonances.
  pyinre(cmn);
  //C
  //C...Reset variables for cross-section calculation.
  FEM_DO_SAFE(i, 0, 200) {
    FEM_DO_SAFE(j, 1, 3) {
      ngen(i, j) = 0;
      xsec(i, j) = 0.f;
    }
  }
  vint(108) = 0.f;
  //C
  //C...Find parametrized total cross-sections.
  if (mint(43) == 4) {
    pyxtot(cmn);
  }
  //C
  //C...Maxima of differential cross-sections.
  if (mstp(121) <= 0) {
    pymaxi(cmn);
  }
  //C
  //C...Initialize possibility of overlayed events.
  if (mstp(131) != 0) {
    pyovly(cmn, 1);
  }
  //C
  //C...Initialize multiple interactions with variable impact parameter.
  if (mint(43) == 4 && (mint(45) != 0 || mstp(131) != 0) && mstp(82) >= 2) {
    pymult(cmn, 1);
  }
  //C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1600)
  //C
  //C...Formats for initialization information.
  //Clin 1000 FORMAT(///20X,'The Lund Monte Carlo - PYTHIA version ',I1,'.',I1/
  //Clin     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)
  //Clin 1100 FORMAT('1',18('*'),1X,'PYINIT: initialization of PYTHIA ',
  //Clin     &'routines',1X,17('*'))
  //Clin 1600 FORMAT(/1X,22('*'),1X,'PYINIT: initialization completed',1X,
  //Clin     &22('*'))
  //C
}

//C
//C*********************************************************************
//C
void
pyspli(
  common& cmn,
  int const& kf,
  int const& kflin,
  int& kflch,
  int& kflsp)
{
  //C
  //C...In case of a hadron remnant which is more complicated than just a
  //C...quark or a diquark, split it into two (partons or hadron + parton).
  //C
  //C...Preliminaries. Parton composition.
  int kfa = fem::iabs(kf);
  int kfs = fem::isign(1, kf);
  arr_1d<3, int> kfl(fem::fill0);
  kfl(1) = fem::mod(kfa / 1000, 10);
  kfl(2) = fem::mod(kfa / 100, 10);
  kfl(3) = fem::mod(kfa / 10, 10);
  int kflr = kflin * kfs;
  kflch = 0;
  //C
  //C...Subdivide meson.
  int kfdump = fem::int0;
  int nagr = fem::int0;
  int j = fem::int0;
  float ragr = fem::float0;
  int iagr = fem::int0;
  int id1 = fem::int0;
  int id2 = fem::int0;
  int ksp = fem::int0;
  if (kfl(1) == 0) {
    kfl(2) = kfl(2) * fem::pow((-1), kfl(2));
    kfl(3) = -kfl(3) * fem::pow((-1), fem::iabs(kfl(2)));
    if (kflr == kfl(2)) {
      kflsp = kfl(3);
    }
    else if (kflr == kfl(3)) {
      kflsp = kfl(2);
    }
    else if (fem::iabs(kflr) == 21 && rlu(cmn, 0) > 0.5f) {
      kflsp = kfl(2);
      kflch = kfl(3);
    }
    else if (fem::iabs(kflr) == 21) {
      kflsp = kfl(3);
      kflch = kfl(2);
    }
    else if (kflr * kfl(2) > 0) {
      lukfdi(cmn, -kflr, kfl(2), kfdump, kflch);
      kflsp = kfl(3);
    }
    else {
      lukfdi(cmn, -kflr, kfl(3), kfdump, kflch);
      kflsp = kfl(2);
    }
    //C
    //C...Subdivide baryon.
  }
  else {
    nagr = 0;
    FEM_DO_SAFE(j, 1, 3) {
      if (kflr == kfl(j)) {
        nagr++;
      }
    }
    if (nagr >= 1) {
      ragr = 0.00001f + (nagr - 0.00002f) * rlu(cmn, 0);
      iagr = 0;
      FEM_DO_SAFE(j, 1, 3) {
        if (kflr == kfl(j)) {
          ragr = ragr - 1.f;
        }
        if (iagr == 0 && ragr <= 0.f) {
          iagr = j;
        }
      }
    }
    else {
      iagr = fem::fint(1.00001f + 2.99998f * rlu(cmn, 0));
    }
    id1 = 1;
    if (iagr == 1) {
      id1 = 2;
    }
    if (iagr == 1 && kfl(3) > kfl(2)) {
      id1 = 3;
    }
    id2 = 6 - iagr - id1;
    ksp = 3;
    if (fem::mod(kfa, 10) == 2 && kfl(1) == kfl(2)) {
      if (iagr != 3 && rlu(cmn, 0) > 0.25f) {
        ksp = 1;
      }
    }
    else if (fem::mod(kfa, 10) == 2 && kfl(2) >= kfl(3)) {
      if (iagr != 1 && rlu(cmn, 0) > 0.25f) {
        ksp = 1;
      }
    }
    else if (fem::mod(kfa, 10) == 2) {
      if (iagr == 1) {
        ksp = 1;
      }
      if (iagr != 1 && rlu(cmn, 0) > 0.75f) {
        ksp = 1;
      }
    }
    kflsp = 1000 * kfl(id1) + 100 * kfl(id2) + ksp;
    if (kflin == 21) {
      kflch = kfl(iagr);
    }
    else if (nagr == 0 && kflr > 0) {
      lukfdi(cmn, -kflr, kfl(iagr), kfdump, kflch);
    }
    else if (nagr == 0) {
      lukfdi(cmn, 10000 + kflsp, -kflr, kfdump, kflch);
      kflsp = kfl(iagr);
    }
  }
  //C
  //C...Add on correct sign for result.
  kflch = kflch * kfs;
  kflsp = kflsp * kfs;
  //C
}

//C
//C*********************************************************************
//C
void
pykcut(
  common& cmn,
  int& mcut)
{
  //C
  //C...Dummy routine, which the user can replace in order to make cuts on
  //C...the kinematics on the parton level before the matrix elements are
  //C...evaluated and the event is generated. The cross-section estimates
  //C...will automatically take these cuts into account, so the given
  //C...values are for the allowed phase space region only. MCUT=0 means
  //C...that the event has passed the cuts, MCUT=1 that it has failed.
  //C
  mcut = 0;
  //C
}

//C
//C*********************************************************************
//C
void
pyrand(
  common& cmn)
{
  common_write write(cmn);
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int> msub(cmn.msub, dimension(200));
  arr_cref<float> ckin(cmn.ckin, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_cref<int, 2> kfpr(cmn.kfpr, dimension(200, 2));
  arr_cref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_cref<int, 2> isig(cmn.isig, dimension(1000, 3));
  arr_cref<float> sigh(cmn.sigh, dimension(1000));
  arr_ref<int, 2> ngen(cmn.ngen, dim1(0, 200).dim2(3));
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  int isub = fem::int0;
  float rsub = fem::float0;
  int i = fem::int0;
  int kfr1 = fem::int0;
  float taur1 = fem::float0;
  float gamr1 = fem::float0;
  int kfr2 = fem::int0;
  float taur2 = fem::float0;
  float gamr2 = fem::float0;
  int is = fem::int0;
  float sh = fem::float0;
  float sqm1 = fem::float0;
  float sqm2 = fem::float0;
  float sqm3 = fem::float0;
  float sqm4 = fem::float0;
  float sqla12 = fem::float0;
  float sqla34 = fem::float0;
  float thter1 = fem::float0;
  float thter2 = fem::float0;
  float thl = fem::float0;
  float thu = fem::float0;
  float thm = fem::float0;
  int jtmax = fem::int0;
  int jt = fem::int0;
  float sqmmin = fem::float0;
  float sqmi = fem::float0;
  float sqmj = fem::float0;
  float sqmf = fem::float0;
  float squa = fem::float0;
  float quar = fem::float0;
  float sqmmax = fem::float0;
  float b = fem::float0;
  float c = fem::float0;
  float expth = fem::float0;
  float tharg = fem::float0;
  float th = fem::float0;
  float ratlog = fem::float0;
  float rtau = fem::float0;
  int mtau = fem::int0;
  float rtaup = fem::float0;
  int mtaup = fem::int0;
  float ryst = fem::float0;
  int myst = fem::int0;
  float rcth = fem::float0;
  int mcth = fem::int0;
  int mcut = fem::int0;
  int nchn = fem::int0;
  float sigs = fem::float0;
  float viol = fem::float0;
  float xdif = fem::float0;
  float rsigs = fem::float0;
  float qt2 = fem::float0;
  float rqqbar = fem::float0;
  int ichn = fem::int0;
  int kfl1 = fem::int0;
  int kfl2 = fem::int0;
  int kfldum = fem::int0;
  //C
  //C...Generates quantities characterizing the high-pT scattering at the
  //C...parton level according to the matrix elements. Chooses incoming,
  //C...reacting partons, their momentum fractions and one of the possible
  //C...subprocesses.
  //C
  //C...Initial values, specifically for (first) semihard interaction.
  mint(17) = 0;
  mint(18) = 0;
  vint(143) = 1.f;
  vint(144) = 1.f;
  if (msub(95) == 1 || mint(82) >= 2) {
    pymult(cmn, 2);
  }
  isub = 0;
  statement_100:
  mint(51) = 0;
  //C
  //C...Choice of process type - first event of overlay.
  if (mint(82) == 1 && (isub <= 90 || isub > 96)) {
    rsub = xsec(0, 1) * rlu(cmn, 0);
    FEM_DO_SAFE(i, 1, 200) {
      if (msub(i) != 1) {
        goto statement_110;
      }
      isub = i;
      rsub = rsub - xsec(i, 1);
      if (rsub <= 0.f) {
        goto statement_120;
      }
      statement_110:;
    }
    statement_120:
    if (isub == 95) {
      isub = 96;
    }
    //C
    //C...Choice of inclusive process type - overlayed events.
  }
  else if (mint(82) >= 2 && isub == 0) {
    rsub = vint(131) * rlu(cmn, 0);
    isub = 96;
    if (rsub > vint(106)) {
      isub = 93;
    }
    if (rsub > vint(106) + vint(104)) {
      isub = 92;
    }
    if (rsub > vint(106) + vint(104) + vint(103)) {
      isub = 91;
    }
  }
  if (mint(82) == 1) {
    ngen(0, 1)++;
  }
  if (mint(82) == 1) {
    ngen(isub, 1)++;
  }
  mint(1) = isub;
  //C
  //C...Find resonances (explicit or implicit in cross-section).
  mint(72) = 0;
  kfr1 = 0;
  if (iset(isub) == 1 || iset(isub) == 3) {
    kfr1 = kfpr(isub, 1);
  }
  else if (isub >= 71 && isub <= 77) {
    kfr1 = 25;
  }
  if (kfr1 != 0) {
    taur1 = fem::pow2(pmas(kfr1, 1)) / vint(2);
    gamr1 = pmas(kfr1, 1) * pmas(kfr1, 2) / vint(2);
    mint(72) = 1;
    mint(73) = kfr1;
    vint(73) = taur1;
    vint(74) = gamr1;
  }
  if (isub == 141) {
    kfr2 = 23;
    taur2 = fem::pow2(pmas(kfr2, 1)) / vint(2);
    gamr2 = pmas(kfr2, 1) * pmas(kfr2, 2) / vint(2);
    mint(72) = 2;
    mint(74) = kfr2;
    vint(75) = taur2;
    vint(76) = gamr2;
  }
  //C
  //C...Find product masses and minimum pT of process,
  //C...optionally with broadening according to a truncated Breit-Wigner.
  vint(63) = 0.f;
  vint(64) = 0.f;
  mint(71) = 0;
  vint(71) = ckin(3);
  if (mint(82) >= 2) {
    vint(71) = 0.f;
  }
  if (iset(isub) == 2 || iset(isub) == 4) {
    FEM_DO_SAFE(i, 1, 2) {
      if (kfpr(isub, i) == 0) {
      }
      else if (mstp(42) <= 0) {
        vint(62 + i) = fem::pow2(pmas(kfpr(isub, i), 1));
      }
      else {
        vint(62 + i) = fem::pow2(ulmass(cmn, kfpr(isub, i)));
      }
    }
    if (fem::min(vint(63), vint(64)) < fem::pow2(ckin(6))) {
      mint(71) = 1;
    }
    if (mint(71) == 1) {
      vint(71) = fem::max(ckin(3), ckin(5));
    }
  }
  //C
  if (iset(isub) == 0) {
    //C...Double or single diffractive, or elastic scattering:
    //C...choose m^2 according to 1/m^2 (diffractive), constant (elastic)
    is = fem::fint(1.5f + rlu(cmn, 0));
    vint(63) = fem::pow2(vint(3));
    vint(64) = fem::pow2(vint(4));
    if (isub == 92 || isub == 93) {
      vint(62 + is) = fem::pow2(parp(111));
    }
    if (isub == 93) {
      vint(65 - is) = fem::pow2(parp(111));
    }
    sh = vint(2);
    sqm1 = fem::pow2(vint(3));
    sqm2 = fem::pow2(vint(4));
    sqm3 = vint(63);
    sqm4 = vint(64);
    sqla12 = fem::pow2((sh - sqm1 - sqm2)) - 4.f * sqm1 * sqm2;
    sqla34 = fem::pow2((sh - sqm3 - sqm4)) - 4.f * sqm3 * sqm4;
    thter1 = sqm1 + sqm2 + sqm3 + sqm4 - (sqm1 - sqm2) * (sqm3 -
      sqm4) / sh - sh;
    thter2 = fem::sqrt(fem::max(0.f, sqla12)) * fem::sqrt(fem::max(0.f,
      sqla34)) / sh;
    thl = 0.5f * (thter1 - thter2);
    thu = 0.5f * (thter1 + thter2);
    thm = fem::min(fem::max(thl, parp(101)), thu);
    jtmax = 0;
    if (isub == 92 || isub == 93) {
      jtmax = isub - 91;
    }
    FEM_DO_SAFE(jt, 1, jtmax) {
      mint(13 + 3 * jt - is * (2 * jt - 3)) = 1;
      sqmmin = vint(59 + 3 * jt - is * (2 * jt - 3));
      sqmi = fem::pow2(vint(8 - 3 * jt + is * (2 * jt - 3)));
      sqmj = fem::pow2(vint(3 * jt - 1 - is * (2 * jt - 3)));
      sqmf = vint(68 - 3 * jt + is * (2 * jt - 3));
      squa = 0.5f * sh / sqmi * ((1.f + (sqmi - sqmj) / sh) * thm +
        sqmi - sqmf - fem::pow2(sqmj) / sh + (sqmi + sqmj) * sqmf /
        sh + fem::pow2((sqmi - sqmj)) / fem::pow2(sh) * sqmf);
      quar = sh / sqmi * (thm * (thm + sh - sqmi - sqmj - sqmf * (
        1.f - (sqmi - sqmj) / sh)) + sqmi * sqmj - sqmj * sqmf * (
        1.f + (sqmi - sqmj - sqmf) / sh));
      sqmmax = squa + fem::sqrt(fem::max(0.f, fem::pow2(squa) - quar));
      if (fem::abs(quar / fem::pow2(squa)) < 1.e-06f) {
        sqmmax = 0.5f * quar / squa;
      }
      sqmmax = fem::min(sqmmax, fem::pow2((vint(1) - fem::sqrt(sqmf))));
      vint(59 + 3 * jt - is * (2 * jt - 3)) = sqmmin * fem::pow((
        sqmmax / sqmmin), rlu(cmn, 0));
    }
    //C...Choose t-hat according to exp(B*t-hat+C*t-hat^2).
    sqm3 = vint(63);
    sqm4 = vint(64);
    sqla34 = fem::pow2((sh - sqm3 - sqm4)) - 4.f * sqm3 * sqm4;
    thter1 = sqm1 + sqm2 + sqm3 + sqm4 - (sqm1 - sqm2) * (sqm3 -
      sqm4) / sh - sh;
    thter2 = fem::sqrt(fem::max(0.f, sqla12)) * fem::sqrt(fem::max(0.f,
      sqla34)) / sh;
    thl = 0.5f * (thter1 - thter2);
    thu = 0.5f * (thter1 + thter2);
    b = vint(121);
    c = vint(122);
    if (isub == 92 || isub == 93) {
      b = 0.5f * b;
      c = 0.5f * c;
    }
    thm = fem::min(fem::max(thl, parp(101)), thu);
    expth = 0.f;
    tharg = b * (thm - thu);
    if (tharg >  - 20.f) {
      expth = fem::exp(tharg);
    }
    statement_150:
    th = thu + fem::log(expth + (1.f - expth) * rlu(cmn, 0)) / b;
    th = fem::max(thm, fem::min(thu, th));
    ratlog = fem::min((b + c * (th + thm)) * (th - thm), (b + c * (
      th + thu)) * (th - thu));
    if (ratlog < fem::log(rlu(cmn, 0))) {
      goto statement_150;
    }
    vint(21) = 1.f;
    vint(22) = 0.f;
    vint(23) = fem::min(1.f, fem::max(-1.f, (2.f * th - thter1) / thter2));
    //C
    //C...Note: in the following, by In is meant the integral over the
    //C...quantity multiplying coefficient cn.
    //C...Choose tau according to h1(tau)/tau, where
    //C...h1(tau) = c0 + I0/I1*c1*1/tau + I0/I2*c2*1/(tau+tau_R) +
    //C...I0/I3*c3*tau/((s*tau-m^2)^2+(m*Gamma)^2) +
    //C...I0/I4*c4*1/(tau+tau_R') +
    //C...I0/I5*c5*tau/((s*tau-m'^2)^2+(m'*Gamma')^2), and
    //C...c0 + c1 + c2 + c3 + c4 + c5 = 1
  }
  else if (iset(isub) >= 1 && iset(isub) <= 4) {
    pyklim(cmn, 1);
    if (mint(51) != 0) {
      goto statement_100;
    }
    rtau = rlu(cmn, 0);
    mtau = 1;
    if (rtau > coef(isub, 1)) {
      mtau = 2;
    }
    if (rtau > coef(isub, 1) + coef(isub, 2)) {
      mtau = 3;
    }
    if (rtau > coef(isub, 1) + coef(isub, 2) + coef(isub, 3)) {
      mtau = 4;
    }
    if (rtau > coef(isub, 1) + coef(isub, 2) + coef(isub, 3) + coef(isub, 4)) {
      mtau = 5;
    }
    if (rtau > coef(isub, 1) + coef(isub, 2) + coef(isub, 3) + coef(isub,
        4) + coef(isub, 5)) {
      mtau = 6;
    }
    pykmap(cmn, 1, mtau, rlu(cmn, 0));
    //C
    //C...2 -> 3, 4 processes:
    //C...Choose tau' according to h4(tau,tau')/tau', where
    //C...h4(tau,tau') = c0 + I0/I1*c1*(1 - tau/tau')^3/tau', and
    //C...c0 + c1 = 1.
    if (iset(isub) == 3 || iset(isub) == 4) {
      pyklim(cmn, 4);
      if (mint(51) != 0) {
        goto statement_100;
      }
      rtaup = rlu(cmn, 0);
      mtaup = 1;
      if (rtaup > coef(isub, 15)) {
        mtaup = 2;
      }
      pykmap(cmn, 4, mtaup, rlu(cmn, 0));
    }
    //C
    //C...Choose y* according to h2(y*), where
    //C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +
    //C...I0/I3*c3*1/cosh(y*), I0 = y*max-y*min, and c1 + c2 + c3 = 1.
    pyklim(cmn, 2);
    if (mint(51) != 0) {
      goto statement_100;
    }
    ryst = rlu(cmn, 0);
    myst = 1;
    if (ryst > coef(isub, 7)) {
      myst = 2;
    }
    if (ryst > coef(isub, 7) + coef(isub, 8)) {
      myst = 3;
    }
    pykmap(cmn, 2, myst, rlu(cmn, 0));
    //C
    //C...2 -> 2 processes:
    //C...Choose cos(theta-hat) (cth) according to h3(cth), where
    //C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +
    //C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,
    //C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products),
    //C...and c0 + c1 + c2 + c3 + c4 = 1.
    pyklim(cmn, 3);
    if (mint(51) != 0) {
      goto statement_100;
    }
    if (iset(isub) == 2 || iset(isub) == 4) {
      rcth = rlu(cmn, 0);
      mcth = 1;
      if (rcth > coef(isub, 10)) {
        mcth = 2;
      }
      if (rcth > coef(isub, 10) + coef(isub, 11)) {
        mcth = 3;
      }
      if (rcth > coef(isub, 10) + coef(isub, 11) + coef(isub, 12)) {
        mcth = 4;
      }
      if (rcth > coef(isub, 10) + coef(isub, 11) + coef(isub, 12) + coef(isub,
          13)) {
        mcth = 5;
      }
      pykmap(cmn, 3, mcth, rlu(cmn, 0));
    }
    //C
    //C...Low-pT or multiple interactions (first semihard interaction).
  }
  else if (iset(isub) == 5) {
    pymult(cmn, 3);
    isub = mint(1);
  }
  //C
  //C...Choose azimuthal angle.
  vint(24) = paru(2) * rlu(cmn, 0);
  //C
  //C...Check against user cuts on kinematics at parton level.
  mint(51) = 0;
  if (isub <= 90 || isub > 100) {
    pyklim(cmn, 0);
  }
  if (mint(51) != 0) {
    goto statement_100;
  }
  if (mint(82) == 1 && mstp(141) >= 1) {
    mcut = 0;
    if (msub(91) + msub(92) + msub(93) + msub(94) + msub(95) == 0) {
      pykcut(cmn, mcut);
    }
    if (mcut != 0) {
      goto statement_100;
    }
  }
  //C
  //C...Calculate differential cross-section for different subprocesses.
  pysigh(cmn, nchn, sigs);
  //C
  //C...Calculations for Monte Carlo estimate of all cross-sections.
  if (mint(82) == 1 && isub <= 90 || isub >= 96) {
    xsec(isub, 2) += sigs;
  }
  else if (mint(82) == 1) {
    xsec(isub, 2) += xsec(isub, 1);
  }
  //C
  //C...Multiple interactions: store results of cross-section calculation.
  if (mint(43) == 4 && mstp(82) >= 3) {
    vint(153) = sigs;
    pymult(cmn, 4);
  }
  //C
  //C...Weighting using estimate of maximum of differential cross-section.
  viol = sigs / xsec(isub, 1);
  if (viol < rlu(cmn, 0)) {
    goto statement_100;
  }
  //C
  //C...Check for possible violation of estimated maximum of differential
  //C...cross-section used in weighting.
  if (mstp(123) <= 0) {
    if (viol > 1.f) {
      write(mstu(11),
        "(1x,'Error: maximum violated by',1p,e11.3,1x,'in event',1x,i7,'.',/,"
        "1x,'Execution stopped!')"),
        viol, ngen(0, 3) + 1;
      write(mstu(11),
        "(1x,'ISUB = ',i3,'; Point of violation:',/,1x,'tau=',1p,e11.3,"
        "', y* =',e11.3,', cthe = ',0p,f11.7,', tau'' =',1p,e11.3)"),
        isub, vint(21), vint(22), vint(23), vint(26);
      FEM_STOP(0);
    }
  }
  else if (mstp(123) == 1) {
    if (viol > vint(108)) {
      vint(108) = viol;
      //C          IF(VIOL.GT.1.) THEN
      //C            WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1
      //C            WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),
      //C     &      VINT(26)
      //C          ENDIF
    }
  }
  else if (viol > vint(108)) {
    vint(108) = viol;
    if (viol > 1.f) {
      xdif = xsec(isub, 1) * (viol - 1.f);
      xsec(isub, 1) += xdif;
      if (msub(isub) == 1 && (isub <= 90 || isub > 96)) {
        xsec(0, 1) += xdif;
      }
      //C          WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1
      //C          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26)
      //C          IF(ISUB.LE.9) THEN
      //C            WRITE(MSTU(11),1300) ISUB,XSEC(ISUB,1)
      //C          ELSEIF(ISUB.LE.99) THEN
      //C            WRITE(MSTU(11),1400) ISUB,XSEC(ISUB,1)
      //C          ELSE
      //C            WRITE(MSTU(11),1500) ISUB,XSEC(ISUB,1)
      //C          ENDIF
      vint(108) = 1.f;
    }
  }
  //C
  //C...Multiple interactions: choose impact parameter.
  vint(148) = 1.f;
  if (mint(43) == 4 && (isub <= 90 || isub >= 96) && mstp(82) >= 3) {
    pymult(cmn, 5);
    if (vint(150) < rlu(cmn, 0)) {
      goto statement_100;
    }
  }
  if (mint(82) == 1 && msub(95) == 1) {
    if (isub <= 90 || isub >= 95) {
      ngen(95, 1)++;
    }
    if (isub <= 90 || isub >= 96) {
      ngen(96, 2)++;
    }
  }
  if (isub <= 90 || isub >= 96) {
    mint(31)++;
  }
  //C
  //C...Choose flavour of reacting partons (and subprocess).
  rsigs = sigs * rlu(cmn, 0);
  qt2 = vint(48);
  rqqbar = parp(87) * (1.f - fem::pow2((qt2 / (qt2 + fem::pow2((parp(
    88) * parp(82)))))));
  if (isub != 95 && (isub != 96 || mstp(82) <= 1 || rlu(cmn, 0) > rqqbar)) {
    FEM_DO_SAFE(ichn, 1, nchn) {
      kfl1 = isig(ichn, 1);
      kfl2 = isig(ichn, 2);
      mint(2) = isig(ichn, 3);
      rsigs = rsigs - sigh(ichn);
      if (rsigs <= 0.f) {
        goto statement_210;
      }
    }
    //C
    //C...Multiple interactions: choose qqbar preferentially at small pT.
  }
  else if (isub == 96) {
    pyspli(cmn, mint(11), 21, kfl1, kfldum);
    pyspli(cmn, mint(12), 21, kfl2, kfldum);
    mint(1) = 11;
    mint(2) = 1;
    if (kfl1 == kfl2 && rlu(cmn, 0) < 0.5f) {
      mint(2) = 2;
    }
    //C
    //C...Low-pT: choose string drawing configuration.
  }
  else {
    kfl1 = 21;
    kfl2 = 21;
    rsigs = 6.f * rlu(cmn, 0);
    mint(2) = 1;
    if (rsigs > 1.f) {
      mint(2) = 2;
    }
    if (rsigs > 2.f) {
      mint(2) = 3;
    }
  }
  //C
  //C...Reassign QCD process. Partons before initial state radiation.
  statement_210:
  if (mint(2) > 10) {
    mint(1) = mint(2) / 10;
    mint(2) = fem::mod(mint(2), 10);
  }
  mint(15) = kfl1;
  mint(16) = kfl2;
  mint(13) = mint(15);
  mint(14) = mint(16);
  vint(141) = vint(41);
  vint(142) = vint(42);
  //C
  //C...Format statements for differential cross-section maximum violations.
  //Clin 1200 FORMAT(1X,'Warning: maximum violated by',1P,E11.3,1X,
  //C     &'in event',1X,I7)
  //C 1300 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,E11.3)
  //C 1400 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,E11.3)
  //Clin 1500 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,E11.3)
  //C
}

//C
//C*********************************************************************
//C
void
pyscat(
  common& cmn)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<float, 2> vckm(cmn.vckm, dimension(4, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_cref<int, 2> kfpr(cmn.kfpr, dimension(200, 2));
  arr_cref<int, 3> icol(cmn.icol, dimension(40, 4, 2));
  //
  int isub = fem::int0;
  int idoc = fem::int0;
  int ipu1 = fem::int0;
  int ipu2 = fem::int0;
  int ipu3 = fem::int0;
  int ipu4 = fem::int0;
  int ipu5 = fem::int0;
  int ipu6 = fem::int0;
  int jt = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  int kfres = fem::int0;
  float sh = fem::float0;
  float shr = fem::float0;
  float shp = fem::float0;
  float shpr = fem::float0;
  float shuser = fem::float0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  arr_1d<41, float> wdtp(dim1(0, 40), fem::fill0);
  arr_2d<41, 6, float> wdte(dim1(0, 40).dim2(0, 5), fem::fill0);
  float rkfl = fem::float0;
  int kflq = fem::int0;
  int js = fem::int0;
  int kcc = fem::int0;
  int kcs = fem::int0;
  int kch1 = fem::int0;
  int kch2 = fem::int0;
  float xh = fem::float0;
  arr_1d<2, float> pmq(fem::fill0);
  float zmin = fem::float0;
  float zmax = fem::float0;
  arr_1d<2, float> z(fem::fill0);
  float sqc1 = fem::float0;
  float c1 = fem::float0;
  float c2 = fem::float0;
  arr_1d<2, float> cthe(fem::fill0);
  float phir = fem::float0;
  float cphi = fem::float0;
  float ang = fem::float0;
  float z1 = fem::float0;
  float z2 = fem::float0;
  float z3 = fem::float0;
  int ia = fem::int0;
  float rvckm = fem::float0;
  int ib = fem::int0;
  int ipm = fem::int0;
  int idc = fem::int0;
  int ja = fem::int0;
  int kfa1 = fem::int0;
  int kfa2 = fem::int0;
  arr_1d<2, float> phi(fem::fill0);
  float pabs = fem::float0;
  float ptabs = fem::float0;
  int izw = fem::int0;
  int ipu = fem::int0;
  float bexcm = fem::float0;
  float beycm = fem::float0;
  float bezcm = fem::float0;
  float gamcm = fem::float0;
  float bepcm = fem::float0;
  float px = fem::float0;
  float py = fem::float0;
  float pz = fem::float0;
  float thecm = fem::float0;
  float phicm = fem::float0;
  float sqlam = fem::float0;
  float cthwz = fem::float0;
  float sthwz = fem::float0;
  float phiwz = fem::float0;
  int jc = fem::int0;
  //C
  //C...Finds outgoing flavours and event type; sets up the kinematics
  //C...and colour flow of the hard scattering.
  //C
  //C...Choice of subprocess, number of documentation lines.
  isub = mint(1);
  idoc = 6 + iset(isub);
  if (isub == 95) {
    idoc = 8;
  }
  mint(3) = idoc - 6;
  if (idoc >= 9) {
    idoc += 2;
  }
  mint(4) = idoc;
  ipu1 = mint(84) + 1;
  ipu2 = mint(84) + 2;
  ipu3 = mint(84) + 3;
  ipu4 = mint(84) + 4;
  ipu5 = mint(84) + 5;
  ipu6 = mint(84) + 6;
  //C
  //C...Reset K, P and V vectors. Store incoming particles.
  FEM_DO_SAFE(jt, 1, mstp(126) + 10) {
    i = mint(83) + jt;
    FEM_DO_SAFE(j, 1, 5) {
      k(i, j) = 0;
      p(i, j) = 0.f;
      v(i, j) = 0.f;
    }
  }
  FEM_DO_SAFE(jt, 1, 2) {
    i = mint(83) + jt;
    k(i, 1) = 21;
    k(i, 2) = mint(10 + jt);
    p(i, 1) = 0.f;
    p(i, 2) = 0.f;
    p(i, 5) = vint(2 + jt);
    p(i, 3) = vint(5) * fem::pow((-1), (jt + 1));
    p(i, 4) = fem::sqrt(fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
  }
  mint(6) = 2;
  kfres = 0;
  //C
  //C...Store incoming partons in their CM-frame.
  sh = vint(44);
  shr = fem::sqrt(sh);
  shp = vint(26) * vint(2);
  shpr = fem::sqrt(shp);
  shuser = shr;
  if (iset(isub) >= 3) {
    shuser = shpr;
  }
  FEM_DO_SAFE(jt, 1, 2) {
    i = mint(84) + jt;
    k(i, 1) = 14;
    k(i, 2) = mint(14 + jt);
    k(i, 3) = mint(83) + 2 + jt;
    p(i, 5) = ulmass(cmn, k(i, 2));
  }
  if (p(ipu1, 5) + p(ipu2, 5) >= shuser) {
    p(ipu1, 5) = 0.f;
    p(ipu2, 5) = 0.f;
  }
  p(ipu1, 4) = 0.5f * (shuser + (fem::pow2(p(ipu1, 5)) - fem::pow2(p(ipu2,
    5))) / shuser);
  p(ipu1, 3) = fem::sqrt(fem::max(0.f, fem::pow2(p(ipu1, 4)) -
    fem::pow2(p(ipu1, 5))));
  p(ipu2, 4) = shuser - p(ipu1, 4);
  p(ipu2, 3) = -p(ipu1, 3);
  //C
  //C...Copy incoming partons to documentation lines.
  FEM_DO_SAFE(jt, 1, 2) {
    i1 = mint(83) + 4 + jt;
    i2 = mint(84) + jt;
    k(i1, 1) = 21;
    k(i1, 2) = k(i2, 2);
    k(i1, 3) = i1 - 2;
    FEM_DO_SAFE(j, 1, 5) {
      p(i1, j) = p(i2, j);
    }
  }
  //C
  //C...Choose new quark flavour for relevant annihilation graphs.
  if (isub == 12 || isub == 53) {
    pywidt(cmn, 21, shr, wdtp, wdte);
    rkfl = (wdte(0, 1) + wdte(0, 2) + wdte(0, 4)) * rlu(cmn, 0);
    FEM_DO_SAFE(i, 1, 2 * mstp(1)) {
      kflq = i;
      rkfl = rkfl - (wdte(i, 1) + wdte(i, 2) + wdte(i, 4));
      if (rkfl <= 0.f) {
        goto statement_150;
      }
    }
    statement_150:;
  }
  //C
  //C...Final state flavours and colour flow: default values.
  js = 1;
  mint(21) = mint(15);
  mint(22) = mint(16);
  mint(23) = 0;
  mint(24) = 0;
  kcc = 20;
  kcs = fem::isign(1, mint(15));
  //C
  if (isub <= 10) {
    if (isub == 1) {
      //C...f + fb -> gamma*/Z0.
      kfres = 23;
      //C
    }
    else if (isub == 2) {
      //C...f + fb' -> W+/- .
      kch1 = kchg(fem::iabs(mint(15)), 1) * fem::isign(1, mint(15));
      kch2 = kchg(fem::iabs(mint(16)), 1) * fem::isign(1, mint(16));
      kfres = fem::isign(24, kch1 + kch2);
      //C
    }
    else if (isub == 3) {
      //C...f + fb -> H0.
      kfres = 25;
      //C
    }
    else if (isub == 4) {
      //C...gamma + W+/- -> W+/-.
      //C
    }
    else if (isub == 5) {
      //C...Z0 + Z0 -> H0.
      xh = sh / shp;
      mint(21) = mint(15);
      mint(22) = mint(16);
      pmq(1) = ulmass(cmn, mint(21));
      pmq(2) = ulmass(cmn, mint(22));
      statement_240:
      jt = fem::fint(1.5f + rlu(cmn, 0));
      zmin = 2.f * pmq(jt) / shpr;
      zmax = 1.f - pmq(3 - jt) / shpr - (sh - fem::pow2(pmq(jt))) / (
        shpr * (shpr - pmq(3 - jt)));
      zmax = fem::min(1.f - xh, zmax);
      z(jt) = zmin + (zmax - zmin) * rlu(cmn, 0);
      if (-1.f + (1.f + xh) / (1.f - z(jt)) - xh / fem::pow2((1.f - z(
          jt))) < fem::pow2((1.f - xh)) / (4.f * xh) * rlu(cmn, 0)) {
        goto statement_240;
      }
      sqc1 = 1.f - 4.f * fem::pow2(pmq(jt)) / (fem::pow2(z(jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_240;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(23, 1)) - fem::pow2(pmq(
        jt))) / (z(jt) * shp);
      cthe(jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (2.f * rlu(cmn,
        0) - 1.f) * c1)) / c1;
      cthe(jt) = fem::min(1.f, fem::max(-1.f, cthe(jt)));
      z(3 - jt) = 1.f - xh / (1.f - z(jt));
      sqc1 = 1.f - 4.f * fem::pow2(pmq(3 - jt)) / (fem::pow2(z(3 - jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_240;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(23, 1)) - fem::pow2(pmq(3 -
        jt))) / (z(3 - jt) * shp);
      cthe(3 - jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (
        2.f * rlu(cmn, 0) - 1.f) * c1)) / c1;
      cthe(3 - jt) = fem::min(1.f, fem::max(-1.f, cthe(3 - jt)));
      phir = paru(2) * rlu(cmn, 0);
      cphi = fem::cos(phir);
      ang = cthe(1) * cthe(2) - fem::sqrt(1.f - fem::pow2(cthe(1))) *
        fem::sqrt(1.f - fem::pow2(cthe(2))) * cphi;
      z1 = 2.f - z(jt);
      z2 = ang * fem::sqrt(fem::pow2(z(jt)) - 4.f * fem::pow2(pmq(jt)) / shp);
      z3 = 1.f - z(jt) - xh + (fem::pow2(pmq(1)) + fem::pow2(pmq(2))) / shp;
      z(3 - jt) = 2.f / (fem::pow2(z1) - fem::pow2(z2)) * (z1 * z3 +
        z2 * fem::sqrt(fem::pow2(z3) - (fem::pow2(z1) - fem::pow2(
        z2)) * fem::pow2(pmq(3 - jt)) / shp));
      zmin = 2.f * pmq(3 - jt) / shpr;
      zmax = 1.f - pmq(jt) / shpr - (sh - fem::pow2(pmq(3 - jt))) / (
        shpr * (shpr - pmq(jt)));
      zmax = fem::min(1.f - xh, zmax);
      if (z(3 - jt) < zmin || z(3 - jt) > zmax) {
        goto statement_240;
      }
      kcc = 22;
      kfres = 25;
      //C
    }
    else if (isub == 6) {
      //C...Z0 + W+/- -> W+/-.
      //C
    }
    else if (isub == 7) {
      //C...W+ + W- -> Z0.
      //C
    }
    else if (isub == 8) {
      //C...W+ + W- -> H0.
      xh = sh / shp;
      statement_250:
      FEM_DO_SAFE(jt, 1, 2) {
        i = mint(14 + jt);
        ia = fem::iabs(i);
        if (ia <= 10) {
          rvckm = vint(180 + i) * rlu(cmn, 0);
          FEM_DO_SAFE(j, 1, mstp(1)) {
            ib = 2 * j - 1 + fem::mod(ia, 2);
            ipm = (5 - fem::isign(1, i)) / 2;
            idc = j + mdcy(ia, 2) + 2;
            if (mdme(idc, 1) != 1 && mdme(idc, 1) != ipm) {
              goto statement_270;
            }
            mint(20 + jt) = fem::isign(ib, i);
            rvckm = rvckm - vckm((ia + 1) / 2, (ib + 1) / 2);
            if (rvckm <= 0.f) {
              goto statement_280;
            }
            statement_270:;
          }
        }
        else {
          ib = 2 * ((ia + 1) / 2) - 1 + fem::mod(ia, 2);
          mint(20 + jt) = fem::isign(ib, i);
        }
        statement_280:
        pmq(jt) = ulmass(cmn, mint(20 + jt));
      }
      jt = fem::fint(1.5f + rlu(cmn, 0));
      zmin = 2.f * pmq(jt) / shpr;
      zmax = 1.f - pmq(3 - jt) / shpr - (sh - fem::pow2(pmq(jt))) / (
        shpr * (shpr - pmq(3 - jt)));
      zmax = fem::min(1.f - xh, zmax);
      z(jt) = zmin + (zmax - zmin) * rlu(cmn, 0);
      if (-1.f + (1.f + xh) / (1.f - z(jt)) - xh / fem::pow2((1.f - z(
          jt))) < fem::pow2((1.f - xh)) / (4.f * xh) * rlu(cmn, 0)) {
        goto statement_250;
      }
      sqc1 = 1.f - 4.f * fem::pow2(pmq(jt)) / (fem::pow2(z(jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_250;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(24, 1)) - fem::pow2(pmq(
        jt))) / (z(jt) * shp);
      cthe(jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (2.f * rlu(cmn,
        0) - 1.f) * c1)) / c1;
      cthe(jt) = fem::min(1.f, fem::max(-1.f, cthe(jt)));
      z(3 - jt) = 1.f - xh / (1.f - z(jt));
      sqc1 = 1.f - 4.f * fem::pow2(pmq(3 - jt)) / (fem::pow2(z(3 - jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_250;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(24, 1)) - fem::pow2(pmq(3 -
        jt))) / (z(3 - jt) * shp);
      cthe(3 - jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (
        2.f * rlu(cmn, 0) - 1.f) * c1)) / c1;
      cthe(3 - jt) = fem::min(1.f, fem::max(-1.f, cthe(3 - jt)));
      phir = paru(2) * rlu(cmn, 0);
      cphi = fem::cos(phir);
      ang = cthe(1) * cthe(2) - fem::sqrt(1.f - fem::pow2(cthe(1))) *
        fem::sqrt(1.f - fem::pow2(cthe(2))) * cphi;
      z1 = 2.f - z(jt);
      z2 = ang * fem::sqrt(fem::pow2(z(jt)) - 4.f * fem::pow2(pmq(jt)) / shp);
      z3 = 1.f - z(jt) - xh + (fem::pow2(pmq(1)) + fem::pow2(pmq(2))) / shp;
      z(3 - jt) = 2.f / (fem::pow2(z1) - fem::pow2(z2)) * (z1 * z3 +
        z2 * fem::sqrt(fem::pow2(z3) - (fem::pow2(z1) - fem::pow2(
        z2)) * fem::pow2(pmq(3 - jt)) / shp));
      zmin = 2.f * pmq(3 - jt) / shpr;
      zmax = 1.f - pmq(jt) / shpr - (sh - fem::pow2(pmq(3 - jt))) / (
        shpr * (shpr - pmq(jt)));
      zmax = fem::min(1.f - xh, zmax);
      if (z(3 - jt) < zmin || z(3 - jt) > zmax) {
        goto statement_250;
      }
      kcc = 22;
      kfres = 25;
    }
    //C
  }
  else if (isub <= 20) {
    if (isub == 11) {
      //C...f + f' -> f + f'; th = (p(f)-p(f))**2.
      kcc = mint(2);
      if (mint(15) * mint(16) < 0) {
        kcc += 2;
      }
      //C
    }
    else if (isub == 12) {
      //C...f + fb -> f' + fb'; th = (p(f)-p(f'))**2.
      mint(21) = fem::isign(kflq, mint(15));
      mint(22) = -mint(21);
      kcc = 4;
      //C
    }
    else if (isub == 13) {
      //C...f + fb -> g + g; th arbitrary.
      mint(21) = 21;
      mint(22) = 21;
      kcc = mint(2) + 4;
      //C
    }
    else if (isub == 14) {
      //C...f + fb -> g + gam; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 21;
      mint(23 - js) = 22;
      kcc = 17 + js;
      //C
    }
    else if (isub == 15) {
      //C...f + fb -> g + Z0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 21;
      mint(23 - js) = 23;
      kcc = 17 + js;
      //C
    }
    else if (isub == 16) {
      //C...f + fb' -> g + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(fem::iabs(mint(15)), 1) * fem::isign(1, mint(15));
      kch2 = kchg(fem::iabs(mint(16)), 1) * fem::isign(1, mint(16));
      if (mint(15) * (kch1 + kch2) < 0) {
        js = 2;
      }
      mint(20 + js) = 21;
      mint(23 - js) = fem::isign(24, kch1 + kch2);
      kcc = 17 + js;
      //C
    }
    else if (isub == 17) {
      //C...f + fb -> g + H0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 21;
      mint(23 - js) = 25;
      kcc = 17 + js;
      //C
    }
    else if (isub == 18) {
      //C...f + fb -> gamma + gamma; th arbitrary.
      mint(21) = 22;
      mint(22) = 22;
      //C
    }
    else if (isub == 19) {
      //C...f + fb -> gamma + Z0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 22;
      mint(23 - js) = 23;
      //C
    }
    else if (isub == 20) {
      //C...f + fb' -> gamma + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(fem::iabs(mint(15)), 1) * fem::isign(1, mint(15));
      kch2 = kchg(fem::iabs(mint(16)), 1) * fem::isign(1, mint(16));
      if (mint(15) * (kch1 + kch2) < 0) {
        js = 2;
      }
      mint(20 + js) = 22;
      mint(23 - js) = fem::isign(24, kch1 + kch2);
    }
    //C
  }
  else if (isub <= 30) {
    if (isub == 21) {
      //C...f + fb -> gamma + H0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 22;
      mint(23 - js) = 25;
      //C
    }
    else if (isub == 22) {
      //C...f + fb -> Z0 + Z0; th arbitrary.
      mint(21) = 23;
      mint(22) = 23;
      //C
    }
    else if (isub == 23) {
      //C...f + fb' -> Z0 + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(fem::iabs(mint(15)), 1) * fem::isign(1, mint(15));
      kch2 = kchg(fem::iabs(mint(16)), 1) * fem::isign(1, mint(16));
      if (mint(15) * (kch1 + kch2) < 0) {
        js = 2;
      }
      mint(20 + js) = 23;
      mint(23 - js) = fem::isign(24, kch1 + kch2);
      //C
    }
    else if (isub == 24) {
      //C...f + fb -> Z0 + H0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 23;
      mint(23 - js) = 25;
      //C
    }
    else if (isub == 25) {
      //C...f + fb -> W+ + W-; th = (p(f)-p(W-))**2.
      mint(21) = -fem::isign(24, mint(15));
      mint(22) = -mint(21);
      //C
    }
    else if (isub == 26) {
      //C...f + fb' -> W+/- + H0; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(fem::iabs(mint(15)), 1) * fem::isign(1, mint(15));
      kch2 = kchg(fem::iabs(mint(16)), 1) * fem::isign(1, mint(16));
      if (mint(15) * (kch1 + kch2) > 0) {
        js = 2;
      }
      mint(20 + js) = fem::isign(24, kch1 + kch2);
      mint(23 - js) = 25;
      //C
    }
    else if (isub == 27) {
      //C...f + fb -> H0 + H0.
      //C
    }
    else if (isub == 28) {
      //C...f + g -> f + g; th = (p(f)-p(f))**2.
      kcc = mint(2) + 6;
      if (mint(15) == 21) {
        kcc += 2;
      }
      if (mint(15) != 21) {
        kcs = fem::isign(1, mint(15));
      }
      if (mint(16) != 21) {
        kcs = fem::isign(1, mint(16));
      }
      //C
    }
    else if (isub == 29) {
      //C...f + g -> f + gamma; th = (p(f)-p(f))**2.
      if (mint(15) == 21) {
        js = 2;
      }
      mint(23 - js) = 22;
      kcc = 15 + js;
      kcs = fem::isign(1, mint(14 + js));
      //C
    }
    else if (isub == 30) {
      //C...f + g -> f + Z0; th = (p(f)-p(f))**2.
      if (mint(15) == 21) {
        js = 2;
      }
      mint(23 - js) = 23;
      kcc = 15 + js;
      kcs = fem::isign(1, mint(14 + js));
    }
    //C
  }
  else if (isub <= 40) {
    if (isub == 31) {
      //C...f + g -> f' + W+/-; th = (p(f)-p(f'))**2; choose flavour f'.
      if (mint(15) == 21) {
        js = 2;
      }
      i = mint(14 + js);
      ia = fem::iabs(i);
      mint(23 - js) = fem::isign(24, kchg(ia, 1) * i);
      rvckm = vint(180 + i) * rlu(cmn, 0);
      FEM_DO_SAFE(j, 1, mstp(1)) {
        ib = 2 * j - 1 + fem::mod(ia, 2);
        ipm = (5 - fem::isign(1, i)) / 2;
        idc = j + mdcy(ia, 2) + 2;
        if (mdme(idc, 1) != 1 && mdme(idc, 1) != ipm) {
          goto statement_220;
        }
        mint(20 + js) = fem::isign(ib, i);
        rvckm = rvckm - vckm((ia + 1) / 2, (ib + 1) / 2);
        if (rvckm <= 0.f) {
          goto statement_230;
        }
        statement_220:;
      }
      statement_230:
      kcc = 15 + js;
      kcs = fem::isign(1, mint(14 + js));
      //C
    }
    else if (isub == 32) {
      //C...f + g -> f + H0; th = (p(f)-p(f))**2.
      if (mint(15) == 21) {
        js = 2;
      }
      mint(23 - js) = 25;
      kcc = 15 + js;
      kcs = fem::isign(1, mint(14 + js));
      //C
    }
    else if (isub == 33) {
      //C...f + gamma -> f + g.
      //C
    }
    else if (isub == 34) {
      //C...f + gamma -> f + gamma.
      //C
    }
    else if (isub == 35) {
      //C...f + gamma -> f + Z0.
      //C
    }
    else if (isub == 36) {
      //C...f + gamma -> f' + W+/-.
      //C
    }
    else if (isub == 37) {
      //C...f + gamma -> f + H0.
      //C
    }
    else if (isub == 38) {
      //C...f + Z0 -> f + g.
      //C
    }
    else if (isub == 39) {
      //C...f + Z0 -> f + gamma.
      //C
    }
    else if (isub == 40) {
      //C...f + Z0 -> f + Z0.
    }
    //C
  }
  else if (isub <= 50) {
    if (isub == 41) {
      //C...f + Z0 -> f' + W+/-.
      //C
    }
    else if (isub == 42) {
      //C...f + Z0 -> f + H0.
      //C
    }
    else if (isub == 43) {
      //C...f + W+/- -> f' + g.
      //C
    }
    else if (isub == 44) {
      //C...f + W+/- -> f' + gamma.
      //C
    }
    else if (isub == 45) {
      //C...f + W+/- -> f' + Z0.
      //C
    }
    else if (isub == 46) {
      //C...f + W+/- -> f' + W+/-.
      //C
    }
    else if (isub == 47) {
      //C...f + W+/- -> f' + H0.
      //C
    }
    else if (isub == 48) {
      //C...f + H0 -> f + g.
      //C
    }
    else if (isub == 49) {
      //C...f + H0 -> f + gamma.
      //C
    }
    else if (isub == 50) {
      //C...f + H0 -> f + Z0.
    }
    //C
  }
  else if (isub <= 60) {
    if (isub == 51) {
      //C...f + H0 -> f' + W+/-.
      //C
    }
    else if (isub == 52) {
      //C...f + H0 -> f + H0.
      //C
    }
    else if (isub == 53) {
      //C...g + g -> f + fb; th arbitrary.
      kcs = fem::pow((-1), fem::fint(1.5f + rlu(cmn, 0)));
      mint(21) = fem::isign(kflq, kcs);
      mint(22) = -mint(21);
      kcc = mint(2) + 10;
      //C
    }
    else if (isub == 54) {
      //C...g + gamma -> f + fb.
      //C
    }
    else if (isub == 55) {
      //C...g + Z0 -> f + fb.
      //C
    }
    else if (isub == 56) {
      //C...g + W+/- -> f + fb'.
      //C
    }
    else if (isub == 57) {
      //C...g + H0 -> f + fb.
      //C
    }
    else if (isub == 58) {
      //C...gamma + gamma -> f + fb.
      //C
    }
    else if (isub == 59) {
      //C...gamma + Z0 -> f + fb.
      //C
    }
    else if (isub == 60) {
      //C...gamma + W+/- -> f + fb'.
    }
    //C
  }
  else if (isub <= 70) {
    if (isub == 61) {
      //C...gamma + H0 -> f + fb.
      //C
    }
    else if (isub == 62) {
      //C...Z0 + Z0 -> f + fb.
      //C
    }
    else if (isub == 63) {
      //C...Z0 + W+/- -> f + fb'.
      //C
    }
    else if (isub == 64) {
      //C...Z0 + H0 -> f + fb.
      //C
    }
    else if (isub == 65) {
      //C...W+ + W- -> f + fb.
      //C
    }
    else if (isub == 66) {
      //C...W+/- + H0 -> f + fb'.
      //C
    }
    else if (isub == 67) {
      //C...H0 + H0 -> f + fb.
      //C
    }
    else if (isub == 68) {
      //C...g + g -> g + g; th arbitrary.
      kcc = mint(2) + 12;
      kcs = fem::pow((-1), fem::fint(1.5f + rlu(cmn, 0)));
      //C
    }
    else if (isub == 69) {
      //C...gamma + gamma -> W+ + W-.
      //C
    }
    else if (isub == 70) {
      //C...gamma + W+/- -> gamma + W+/-
    }
    //C
  }
  else if (isub <= 80) {
    if (isub == 71 || isub == 72) {
      //C...Z0 + Z0 -> Z0 + Z0; Z0 + Z0 -> W+ + W-.
      xh = sh / shp;
      mint(21) = mint(15);
      mint(22) = mint(16);
      pmq(1) = ulmass(cmn, mint(21));
      pmq(2) = ulmass(cmn, mint(22));
      statement_290:
      jt = fem::fint(1.5f + rlu(cmn, 0));
      zmin = 2.f * pmq(jt) / shpr;
      zmax = 1.f - pmq(3 - jt) / shpr - (sh - fem::pow2(pmq(jt))) / (
        shpr * (shpr - pmq(3 - jt)));
      zmax = fem::min(1.f - xh, zmax);
      z(jt) = zmin + (zmax - zmin) * rlu(cmn, 0);
      if (-1.f + (1.f + xh) / (1.f - z(jt)) - xh / fem::pow2((1.f - z(
          jt))) < fem::pow2((1.f - xh)) / (4.f * xh) * rlu(cmn, 0)) {
        goto statement_290;
      }
      sqc1 = 1.f - 4.f * fem::pow2(pmq(jt)) / (fem::pow2(z(jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_290;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(23, 1)) - fem::pow2(pmq(
        jt))) / (z(jt) * shp);
      cthe(jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (2.f * rlu(cmn,
        0) - 1.f) * c1)) / c1;
      cthe(jt) = fem::min(1.f, fem::max(-1.f, cthe(jt)));
      z(3 - jt) = 1.f - xh / (1.f - z(jt));
      sqc1 = 1.f - 4.f * fem::pow2(pmq(3 - jt)) / (fem::pow2(z(3 - jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_290;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(23, 1)) - fem::pow2(pmq(3 -
        jt))) / (z(3 - jt) * shp);
      cthe(3 - jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (
        2.f * rlu(cmn, 0) - 1.f) * c1)) / c1;
      cthe(3 - jt) = fem::min(1.f, fem::max(-1.f, cthe(3 - jt)));
      phir = paru(2) * rlu(cmn, 0);
      cphi = fem::cos(phir);
      ang = cthe(1) * cthe(2) - fem::sqrt(1.f - fem::pow2(cthe(1))) *
        fem::sqrt(1.f - fem::pow2(cthe(2))) * cphi;
      z1 = 2.f - z(jt);
      z2 = ang * fem::sqrt(fem::pow2(z(jt)) - 4.f * fem::pow2(pmq(jt)) / shp);
      z3 = 1.f - z(jt) - xh + (fem::pow2(pmq(1)) + fem::pow2(pmq(2))) / shp;
      z(3 - jt) = 2.f / (fem::pow2(z1) - fem::pow2(z2)) * (z1 * z3 +
        z2 * fem::sqrt(fem::pow2(z3) - (fem::pow2(z1) - fem::pow2(
        z2)) * fem::pow2(pmq(3 - jt)) / shp));
      zmin = 2.f * pmq(3 - jt) / shpr;
      zmax = 1.f - pmq(jt) / shpr - (sh - fem::pow2(pmq(3 - jt))) / (
        shpr * (shpr - pmq(jt)));
      zmax = fem::min(1.f - xh, zmax);
      if (z(3 - jt) < zmin || z(3 - jt) > zmax) {
        goto statement_290;
      }
      kcc = 22;
      //C
    }
    else if (isub == 73) {
      //C...Z0 + W+/- -> Z0 + W+/-.
      xh = sh / shp;
      statement_300:
      jt = fem::fint(1.5f + rlu(cmn, 0));
      i = mint(14 + jt);
      ia = fem::iabs(i);
      if (ia <= 10) {
        rvckm = vint(180 + i) * rlu(cmn, 0);
        FEM_DO_SAFE(j, 1, mstp(1)) {
          ib = 2 * j - 1 + fem::mod(ia, 2);
          ipm = (5 - fem::isign(1, i)) / 2;
          idc = j + mdcy(ia, 2) + 2;
          if (mdme(idc, 1) != 1 && mdme(idc, 1) != ipm) {
            goto statement_320;
          }
          mint(20 + jt) = fem::isign(ib, i);
          rvckm = rvckm - vckm((ia + 1) / 2, (ib + 1) / 2);
          if (rvckm <= 0.f) {
            goto statement_330;
          }
          statement_320:;
        }
      }
      else {
        ib = 2 * ((ia + 1) / 2) - 1 + fem::mod(ia, 2);
        mint(20 + jt) = fem::isign(ib, i);
      }
      statement_330:
      pmq(jt) = ulmass(cmn, mint(20 + jt));
      mint(23 - jt) = mint(17 - jt);
      pmq(3 - jt) = ulmass(cmn, mint(23 - jt));
      jt = fem::fint(1.5f + rlu(cmn, 0));
      zmin = 2.f * pmq(jt) / shpr;
      zmax = 1.f - pmq(3 - jt) / shpr - (sh - fem::pow2(pmq(jt))) / (
        shpr * (shpr - pmq(3 - jt)));
      zmax = fem::min(1.f - xh, zmax);
      z(jt) = zmin + (zmax - zmin) * rlu(cmn, 0);
      if (-1.f + (1.f + xh) / (1.f - z(jt)) - xh / fem::pow2((1.f - z(
          jt))) < fem::pow2((1.f - xh)) / (4.f * xh) * rlu(cmn, 0)) {
        goto statement_300;
      }
      sqc1 = 1.f - 4.f * fem::pow2(pmq(jt)) / (fem::pow2(z(jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_300;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(23, 1)) - fem::pow2(pmq(
        jt))) / (z(jt) * shp);
      cthe(jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (2.f * rlu(cmn,
        0) - 1.f) * c1)) / c1;
      cthe(jt) = fem::min(1.f, fem::max(-1.f, cthe(jt)));
      z(3 - jt) = 1.f - xh / (1.f - z(jt));
      sqc1 = 1.f - 4.f * fem::pow2(pmq(3 - jt)) / (fem::pow2(z(3 - jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_300;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(23, 1)) - fem::pow2(pmq(3 -
        jt))) / (z(3 - jt) * shp);
      cthe(3 - jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (
        2.f * rlu(cmn, 0) - 1.f) * c1)) / c1;
      cthe(3 - jt) = fem::min(1.f, fem::max(-1.f, cthe(3 - jt)));
      phir = paru(2) * rlu(cmn, 0);
      cphi = fem::cos(phir);
      ang = cthe(1) * cthe(2) - fem::sqrt(1.f - fem::pow2(cthe(1))) *
        fem::sqrt(1.f - fem::pow2(cthe(2))) * cphi;
      z1 = 2.f - z(jt);
      z2 = ang * fem::sqrt(fem::pow2(z(jt)) - 4.f * fem::pow2(pmq(jt)) / shp);
      z3 = 1.f - z(jt) - xh + (fem::pow2(pmq(1)) + fem::pow2(pmq(2))) / shp;
      z(3 - jt) = 2.f / (fem::pow2(z1) - fem::pow2(z2)) * (z1 * z3 +
        z2 * fem::sqrt(fem::pow2(z3) - (fem::pow2(z1) - fem::pow2(
        z2)) * fem::pow2(pmq(3 - jt)) / shp));
      zmin = 2.f * pmq(3 - jt) / shpr;
      zmax = 1.f - pmq(jt) / shpr - (sh - fem::pow2(pmq(3 - jt))) / (
        shpr * (shpr - pmq(jt)));
      zmax = fem::min(1.f - xh, zmax);
      if (z(3 - jt) < zmin || z(3 - jt) > zmax) {
        goto statement_300;
      }
      kcc = 22;
      //C
    }
    else if (isub == 74) {
      //C...Z0 + H0 -> Z0 + H0.
      //C
    }
    else if (isub == 75) {
      //C...W+ + W- -> gamma + gamma.
      //C
    }
    else if (isub == 76 || isub == 77) {
      //C...W+ + W- -> Z0 + Z0; W+ + W- -> W+ + W-.
      xh = sh / shp;
      statement_340:
      FEM_DO_SAFE(jt, 1, 2) {
        i = mint(14 + jt);
        ia = fem::iabs(i);
        if (ia <= 10) {
          rvckm = vint(180 + i) * rlu(cmn, 0);
          FEM_DO_SAFE(j, 1, mstp(1)) {
            ib = 2 * j - 1 + fem::mod(ia, 2);
            ipm = (5 - fem::isign(1, i)) / 2;
            idc = j + mdcy(ia, 2) + 2;
            if (mdme(idc, 1) != 1 && mdme(idc, 1) != ipm) {
              goto statement_360;
            }
            mint(20 + jt) = fem::isign(ib, i);
            rvckm = rvckm - vckm((ia + 1) / 2, (ib + 1) / 2);
            if (rvckm <= 0.f) {
              goto statement_370;
            }
            statement_360:;
          }
        }
        else {
          ib = 2 * ((ia + 1) / 2) - 1 + fem::mod(ia, 2);
          mint(20 + jt) = fem::isign(ib, i);
        }
        statement_370:
        pmq(jt) = ulmass(cmn, mint(20 + jt));
      }
      jt = fem::fint(1.5f + rlu(cmn, 0));
      zmin = 2.f * pmq(jt) / shpr;
      zmax = 1.f - pmq(3 - jt) / shpr - (sh - fem::pow2(pmq(jt))) / (
        shpr * (shpr - pmq(3 - jt)));
      zmax = fem::min(1.f - xh, zmax);
      z(jt) = zmin + (zmax - zmin) * rlu(cmn, 0);
      if (-1.f + (1.f + xh) / (1.f - z(jt)) - xh / fem::pow2((1.f - z(
          jt))) < fem::pow2((1.f - xh)) / (4.f * xh) * rlu(cmn, 0)) {
        goto statement_340;
      }
      sqc1 = 1.f - 4.f * fem::pow2(pmq(jt)) / (fem::pow2(z(jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_340;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(24, 1)) - fem::pow2(pmq(
        jt))) / (z(jt) * shp);
      cthe(jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (2.f * rlu(cmn,
        0) - 1.f) * c1)) / c1;
      cthe(jt) = fem::min(1.f, fem::max(-1.f, cthe(jt)));
      z(3 - jt) = 1.f - xh / (1.f - z(jt));
      sqc1 = 1.f - 4.f * fem::pow2(pmq(3 - jt)) / (fem::pow2(z(3 - jt)) * shp);
      if (sqc1 < 1.e-8f) {
        goto statement_340;
      }
      c1 = fem::sqrt(sqc1);
      c2 = 1.f + 2.f * (fem::pow2(pmas(24, 1)) - fem::pow2(pmq(3 -
        jt))) / (z(3 - jt) * shp);
      cthe(3 - jt) = (c2 - (fem::pow2(c2) - fem::pow2(c1)) / (c2 + (
        2.f * rlu(cmn, 0) - 1.f) * c1)) / c1;
      cthe(3 - jt) = fem::min(1.f, fem::max(-1.f, cthe(3 - jt)));
      phir = paru(2) * rlu(cmn, 0);
      cphi = fem::cos(phir);
      ang = cthe(1) * cthe(2) - fem::sqrt(1.f - fem::pow2(cthe(1))) *
        fem::sqrt(1.f - fem::pow2(cthe(2))) * cphi;
      z1 = 2.f - z(jt);
      z2 = ang * fem::sqrt(fem::pow2(z(jt)) - 4.f * fem::pow2(pmq(jt)) / shp);
      z3 = 1.f - z(jt) - xh + (fem::pow2(pmq(1)) + fem::pow2(pmq(2))) / shp;
      z(3 - jt) = 2.f / (fem::pow2(z1) - fem::pow2(z2)) * (z1 * z3 +
        z2 * fem::sqrt(fem::pow2(z3) - (fem::pow2(z1) - fem::pow2(
        z2)) * fem::pow2(pmq(3 - jt)) / shp));
      zmin = 2.f * pmq(3 - jt) / shpr;
      zmax = 1.f - pmq(jt) / shpr - (sh - fem::pow2(pmq(3 - jt))) / (
        shpr * (shpr - pmq(jt)));
      zmax = fem::min(1.f - xh, zmax);
      if (z(3 - jt) < zmin || z(3 - jt) > zmax) {
        goto statement_340;
      }
      kcc = 22;
      //C
    }
    else if (isub == 78) {
      //C...W+/- + H0 -> W+/- + H0.
      //C
    }
    else if (isub == 79) {
      //C...H0 + H0 -> H0 + H0.
    }
    //C
  }
  else if (isub <= 90) {
    if (isub == 81) {
      //C...q + qb -> Q' + Qb'; th = (p(q)-p(q'))**2.
      mint(21) = fem::isign(mint(46), mint(15));
      mint(22) = -mint(21);
      kcc = 4;
      //C
    }
    else if (isub == 82) {
      //C...g + g -> Q + Qb; th arbitrary.
      kcs = fem::pow((-1), fem::fint(1.5f + rlu(cmn, 0)));
      mint(21) = fem::isign(mint(46), kcs);
      mint(22) = -mint(21);
      kcc = mint(2) + 10;
    }
    //C
  }
  else if (isub <= 100) {
    if (isub == 95) {
      //C...Low-pT ( = energyless g + g -> g + g).
      kcc = mint(2) + 12;
      kcs = fem::pow((-1), fem::fint(1.5f + rlu(cmn, 0)));
      //C
    }
    else if (isub == 96) {
      //C...Multiple interactions (should be reassigned to QCD process).
    }
    //C
  }
  else if (isub <= 110) {
    if (isub == 101) {
      //C...g + g -> gamma*/Z0.
      kcc = 21;
      kfres = 22;
      //C
    }
    else if (isub == 102) {
      //C...g + g -> H0.
      kcc = 21;
      kfres = 25;
    }
    //C
  }
  else if (isub <= 120) {
    if (isub == 111) {
      //C...f + fb -> g + H0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(20 + js) = 21;
      mint(23 - js) = 25;
      kcc = 17 + js;
      //C
    }
    else if (isub == 112) {
      //C...f + g -> f + H0; th = (p(f) - p(f))**2.
      if (mint(15) == 21) {
        js = 2;
      }
      mint(23 - js) = 25;
      kcc = 15 + js;
      kcs = fem::isign(1, mint(14 + js));
      //C
    }
    else if (isub == 113) {
      //C...g + g -> g + H0; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(23 - js) = 25;
      kcc = 22 + js;
      kcs = fem::pow((-1), fem::fint(1.5f + rlu(cmn, 0)));
      //C
    }
    else if (isub == 114) {
      //C...g + g -> gamma + gamma; th arbitrary.
      if (rlu(cmn, 0) > 0.5f) {
        js = 2;
      }
      mint(21) = 22;
      mint(22) = 22;
      kcc = 21;
      //C
    }
    else if (isub == 115) {
      //C...g + g -> gamma + Z0.
      //C
    }
    else if (isub == 116) {
      //C...g + g -> Z0 + Z0.
      //C
    }
    else if (isub == 117) {
      //C...g + g -> W+ + W-.
    }
    //C
  }
  else if (isub <= 140) {
    if (isub == 121) {
      //C...g + g -> f + fb + H0.
    }
    //C
  }
  else if (isub <= 160) {
    if (isub == 141) {
      //C...f + fb -> gamma*/Z0/Z'0.
      kfres = 32;
      //C
    }
    else if (isub == 142) {
      //C...f + fb' -> H+/-.
      kch1 = kchg(fem::iabs(mint(15)), 1) * fem::isign(1, mint(15));
      kch2 = kchg(fem::iabs(mint(16)), 1) * fem::isign(1, mint(16));
      kfres = fem::isign(37, kch1 + kch2);
      //C
    }
    else if (isub == 143) {
      //C...f + fb' -> R.
      kfres = fem::isign(40, mint(15) + mint(16));
    }
    //C
  }
  else {
    if (isub == 161) {
      //C...g + f -> H+/- + f'; th = (p(f)-p(f))**2.
      if (mint(16) == 21) {
        js = 2;
      }
      ia = fem::iabs(mint(17 - js));
      mint(20 + js) = fem::isign(37, kchg(ia, 1) * mint(17 - js));
      ja = ia + fem::mod(ia, 2) - fem::mod(ia + 1, 2);
      mint(23 - js) = fem::isign(ja, mint(17 - js));
      kcc = 18 - js;
      if (mint(15) != 21) {
        kcs = fem::isign(1, mint(15));
      }
      if (mint(16) != 21) {
        kcs = fem::isign(1, mint(16));
      }
    }
  }
  //C
  if (idoc == 7) {
    //C...Resonance not decaying: store colour connection indices.
    i = mint(83) + 7;
    k(ipu3, 1) = 1;
    k(ipu3, 2) = kfres;
    k(ipu3, 3) = i;
    p(ipu3, 4) = shuser;
    p(ipu3, 5) = shuser;
    k(ipu1, 4) = ipu2;
    k(ipu1, 5) = ipu2;
    k(ipu2, 4) = ipu1;
    k(ipu2, 5) = ipu1;
    k(i, 1) = 21;
    k(i, 2) = kfres;
    p(i, 4) = shuser;
    p(i, 5) = shuser;
    n = ipu3;
    mint(21) = kfres;
    mint(22) = 0;
    //C
  }
  else if (idoc == 8) {
    //C...2 -> 2 processes: store outgoing partons in their CM-frame.
    FEM_DO_SAFE(jt, 1, 2) {
      i = mint(84) + 2 + jt;
      k(i, 1) = 1;
      if (fem::iabs(mint(20 + jt)) <= 10 || mint(20 + jt) == 21) {
        k(i, 1) = 3;
      }
      k(i, 2) = mint(20 + jt);
      k(i, 3) = mint(83) + idoc + jt - 2;
      if (fem::iabs(k(i, 2)) <= 10 || k(i, 2) == 21) {
        p(i, 5) = ulmass(cmn, k(i, 2));
      }
      else {
        p(i, 5) = fem::sqrt(vint(63 + fem::mod(js + jt, 2)));
      }
    }
    if (p(ipu3, 5) + p(ipu4, 5) >= shr) {
      kfa1 = fem::iabs(mint(21));
      kfa2 = fem::iabs(mint(22));
      if ((kfa1 > 3 && kfa1 != 21) || (kfa2 > 3 && kfa2 != 21)) {
        mint(51) = 1;
        return;
      }
      p(ipu3, 5) = 0.f;
      p(ipu4, 5) = 0.f;
    }
    p(ipu3, 4) = 0.5f * (shr + (fem::pow2(p(ipu3, 5)) - fem::pow2(p(ipu4,
      5))) / shr);
    p(ipu3, 3) = fem::sqrt(fem::max(0.f, fem::pow2(p(ipu3, 4)) -
      fem::pow2(p(ipu3, 5))));
    p(ipu4, 4) = shr - p(ipu3, 4);
    p(ipu4, 3) = -p(ipu3, 3);
    n = ipu4;
    mint(7) = mint(83) + 7;
    mint(8) = mint(83) + 8;
    //C
    //C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4).
    ludbrb(ipu3, ipu4, fem::acos(vint(23)), vint(24), 0e0, 0e0, 0e0);
    //C
  }
  else if (idoc == 9) {
    //C'''2 -> 3 processes:
    //C
  }
  else if (idoc == 11) {
    //C...Z0 + Z0 -> H0, W+ + W- -> H0: store Higgs and outgoing partons.
    phi(1) = paru(2) * rlu(cmn, 0);
    phi(2) = phi(1) - phir;
    FEM_DO_SAFE(jt, 1, 2) {
      i = mint(84) + 2 + jt;
      k(i, 1) = 1;
      if (fem::iabs(mint(20 + jt)) <= 10 || mint(20 + jt) == 21) {
        k(i, 1) = 3;
      }
      k(i, 2) = mint(20 + jt);
      k(i, 3) = mint(83) + idoc + jt - 2;
      p(i, 5) = ulmass(cmn, k(i, 2));
      if (0.5f * shpr * z(jt) <= p(i, 5)) {
        p(i, 5) = 0.f;
      }
      pabs = fem::sqrt(fem::max(0.f, fem::pow2((0.5f * shpr * z(
        jt))) - fem::pow2(p(i, 5))));
      ptabs = pabs * fem::sqrt(fem::max(0.f, 1.f - fem::pow2(cthe(jt))));
      p(i, 1) = ptabs * fem::cos(phi(jt));
      p(i, 2) = ptabs * fem::sin(phi(jt));
      p(i, 3) = pabs * cthe(jt) * fem::pow((-1), (jt + 1));
      p(i, 4) = 0.5f * shpr * z(jt);
      izw = mint(83) + 6 + jt;
      k(izw, 1) = 21;
      k(izw, 2) = 23;
      if (isub == 8) {
        k(izw, 2) = fem::isign(24, luchge(cmn, mint(14 + jt)));
      }
      k(izw, 3) = izw - 2;
      p(izw, 1) = -p(i, 1);
      p(izw, 2) = -p(i, 2);
      p(izw, 3) = (0.5f * shpr - pabs * cthe(jt)) * fem::pow((-1), (jt + 1));
      p(izw, 4) = 0.5f * shpr * (1.f - z(jt));
      p(izw, 5) = -fem::sqrt(fem::max(0.f, fem::pow2(p(izw, 3)) +
        fem::pow2(ptabs) - fem::pow2(p(izw, 4))));
    }
    i = mint(83) + 9;
    k(ipu5, 1) = 1;
    k(ipu5, 2) = kfres;
    k(ipu5, 3) = i;
    p(ipu5, 5) = shr;
    p(ipu5, 1) = -p(ipu3, 1) - p(ipu4, 1);
    p(ipu5, 2) = -p(ipu3, 2) - p(ipu4, 2);
    p(ipu5, 3) = -p(ipu3, 3) - p(ipu4, 3);
    p(ipu5, 4) = shpr - p(ipu3, 4) - p(ipu4, 4);
    k(i, 1) = 21;
    k(i, 2) = kfres;
    FEM_DO_SAFE(j, 1, 5) {
      p(i, j) = p(ipu5, j);
    }
    n = ipu5;
    mint(23) = kfres;
    //C
  }
  else if (idoc == 12) {
    //C...Z0 and W+/- scattering: store bosons and outgoing partons.
    phi(1) = paru(2) * rlu(cmn, 0);
    phi(2) = phi(1) - phir;
    FEM_DO_SAFE(jt, 1, 2) {
      i = mint(84) + 2 + jt;
      k(i, 1) = 1;
      if (fem::iabs(mint(20 + jt)) <= 10 || mint(20 + jt) == 21) {
        k(i, 1) = 3;
      }
      k(i, 2) = mint(20 + jt);
      k(i, 3) = mint(83) + idoc + jt - 2;
      p(i, 5) = ulmass(cmn, k(i, 2));
      if (0.5f * shpr * z(jt) <= p(i, 5)) {
        p(i, 5) = 0.f;
      }
      pabs = fem::sqrt(fem::max(0.f, fem::pow2((0.5f * shpr * z(
        jt))) - fem::pow2(p(i, 5))));
      ptabs = pabs * fem::sqrt(fem::max(0.f, 1.f - fem::pow2(cthe(jt))));
      p(i, 1) = ptabs * fem::cos(phi(jt));
      p(i, 2) = ptabs * fem::sin(phi(jt));
      p(i, 3) = pabs * cthe(jt) * fem::pow((-1), (jt + 1));
      p(i, 4) = 0.5f * shpr * z(jt);
      izw = mint(83) + 6 + jt;
      k(izw, 1) = 21;
      if (mint(14 + jt) == mint(20 + jt)) {
        k(izw, 2) = 23;
      }
      else {
        k(izw, 2) = fem::isign(24, luchge(cmn, mint(14 + jt)) - luchge(cmn,
          mint(20 + jt)));
      }
      k(izw, 3) = izw - 2;
      p(izw, 1) = -p(i, 1);
      p(izw, 2) = -p(i, 2);
      p(izw, 3) = (0.5f * shpr - pabs * cthe(jt)) * fem::pow((-1), (jt + 1));
      p(izw, 4) = 0.5f * shpr * (1.f - z(jt));
      p(izw, 5) = -fem::sqrt(fem::max(0.f, fem::pow2(p(izw, 3)) +
        fem::pow2(ptabs) - fem::pow2(p(izw, 4))));
      ipu = mint(84) + 4 + jt;
      k(ipu, 1) = 3;
      k(ipu, 2) = kfpr(isub, jt);
      k(ipu, 3) = mint(83) + 8 + jt;
      if (fem::iabs(k(ipu, 2)) <= 10 || k(ipu, 2) == 21) {
        p(ipu, 5) = ulmass(cmn, k(ipu, 2));
      }
      else {
        p(ipu, 5) = fem::sqrt(vint(63 + fem::mod(js + jt, 2)));
      }
      mint(22 + jt) = k(izw, 2);
    }
    if (isub == 72) {
      k(mint(84) + 4 + fem::fint(1.5f + rlu(cmn, 0)), 2) = -24;
    }
    //C...Find rotation and boost for hard scattering subsystem.
    i1 = mint(83) + 7;
    i2 = mint(83) + 8;
    bexcm = (p(i1, 1) + p(i2, 1)) / (p(i1, 4) + p(i2, 4));
    beycm = (p(i1, 2) + p(i2, 2)) / (p(i1, 4) + p(i2, 4));
    bezcm = (p(i1, 3) + p(i2, 3)) / (p(i1, 4) + p(i2, 4));
    gamcm = (p(i1, 4) + p(i2, 4)) / shr;
    bepcm = bexcm * p(i1, 1) + beycm * p(i1, 2) + bezcm * p(i1, 3);
    px = p(i1, 1) + gamcm * (gamcm / (1.f + gamcm) * bepcm - p(i1, 4)) * bexcm;
    py = p(i1, 2) + gamcm * (gamcm / (1.f + gamcm) * bepcm - p(i1, 4)) * beycm;
    pz = p(i1, 3) + gamcm * (gamcm / (1.f + gamcm) * bepcm - p(i1, 4)) * bezcm;
    thecm = ulangl(cmn, pz, fem::sqrt(fem::pow2(px) + fem::pow2(py)));
    phicm = ulangl(cmn, px, py);
    //C...Store hard scattering subsystem. Rotate and boost it.
    sqlam = fem::pow2((sh - fem::pow2(p(ipu5, 5)) - fem::pow2(p(ipu6,
      5)))) - 4.f * fem::pow2(p(ipu5, 5)) * fem::pow2(p(ipu6, 5));
    pabs = fem::sqrt(fem::max(0.f, sqlam / (4.f * sh)));
    cthwz = vint(23);
    sthwz = fem::sqrt(fem::max(0.f, 1.f - fem::pow2(cthwz)));
    phiwz = vint(24) - phicm;
    p(ipu5, 1) = pabs * sthwz * fem::cos(phiwz);
    p(ipu5, 2) = pabs * sthwz * fem::sin(phiwz);
    p(ipu5, 3) = pabs * cthwz;
    p(ipu5, 4) = fem::sqrt(fem::pow2(pabs) + fem::pow2(p(ipu5, 5)));
    p(ipu6, 1) = -p(ipu5, 1);
    p(ipu6, 2) = -p(ipu5, 2);
    p(ipu6, 3) = -p(ipu5, 3);
    p(ipu6, 4) = fem::sqrt(fem::pow2(pabs) + fem::pow2(p(ipu6, 5)));
    ludbrb(ipu5, ipu6, thecm, phicm, fem::dble(bexcm), fem::dble(beycm),
      fem::dble(bezcm));
    FEM_DO_SAFE(jt, 1, 2) {
      i1 = mint(83) + 8 + jt;
      i2 = mint(84) + 4 + jt;
      k(i1, 1) = 21;
      k(i1, 2) = k(i2, 2);
      FEM_DO_SAFE(j, 1, 5) {
        p(i1, j) = p(i2, j);
      }
    }
    n = ipu6;
    mint(7) = mint(83) + 9;
    mint(8) = mint(83) + 10;
  }
  //C
  if (idoc >= 8) {
    //C...Store colour connection indices.
    FEM_DO_SAFE(j, 1, 2) {
      jc = j;
      if (kcs ==  - 1) {
        jc = 3 - j;
      }
      if (icol(kcc, 1, jc) != 0 && k(ipu1, 1) == 14) {
        k(ipu1, j + 3) += mint(84) + icol(kcc, 1, jc);
      }
      if (icol(kcc, 2, jc) != 0 && k(ipu2, 1) == 14) {
        k(ipu2, j + 3) += mint(84) + icol(kcc, 2, jc);
      }
      if (icol(kcc, 3, jc) != 0 && k(ipu3, 1) == 3) {
        k(ipu3, j + 3) = mstu(5) * (mint(84) + icol(kcc, 3, jc));
      }
      if (icol(kcc, 4, jc) != 0 && k(ipu4, 1) == 3) {
        k(ipu4, j + 3) = mstu(5) * (mint(84) + icol(kcc, 4, jc));
      }
    }
    //C
    //C...Copy outgoing partons to documentation lines.
    FEM_DO_SAFE(i, 1, 2) {
      i1 = mint(83) + idoc - 2 + i;
      i2 = mint(84) + 2 + i;
      k(i1, 1) = 21;
      k(i1, 2) = k(i2, 2);
      if (idoc <= 9) {
        k(i1, 3) = 0;
      }
      if (idoc >= 11) {
        k(i1, 3) = mint(83) + 2 + i;
      }
      FEM_DO_SAFE(j, 1, 5) {
        p(i1, j) = p(i2, j);
      }
    }
  }
  mint(52) = n;
  //C
  //C...Low-pT events: remove gluons used for string drawing purposes.
  if (isub == 95) {
    k(ipu3, 1) += 10;
    k(ipu4, 1) += 10;
    FEM_DO_SAFE(j, 41, 66) {
      vint(j) = 0.f;
    }
    FEM_DO_SAFE(i, mint(83) + 5, mint(83) + 8) {
      FEM_DO_SAFE(j, 1, 5) {
        p(i, j) = 0.f;
      }
    }
  }
  //C
}

//C
//C*********************************************************************
//C
void
pysspa(
  common& cmn,
  int& ipu1,
  int& ipu2)
{
  common_write write(cmn);
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<float> paru(cmn.paru, dimension(200));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_cref<float, 2> xsfx(cmn.xsfx, dim1(2).dim2(-40, 40));
  //
  int ipus1 = fem::int0;
  int ipus2 = fem::int0;
  int isub = fem::int0;
  float q2e = fem::float0;
  float tmax = fem::float0;
  float xe0 = fem::float0;
  float alams = fem::float0;
  int ns = fem::int0;
  int jt = fem::int0;
  arr_1d<4, int> kfls(fem::fill0);
  arr_1d<2, float> xs(fem::fill0);
  arr_1d<2, float> zs(fem::fill0);
  arr_1d<2, float> q2s(fem::fill0);
  arr_1d<2, float> tevs(fem::fill0);
  arr_1d<2, float> alam(fem::fill0);
  arr_1d<2, float> the2(fem::fill0);
  int kfl = fem::int0;
  arr_2d<2, 13, float> xfs(dim1(2).dim2(-6, 6), fem::fill0);
  double dsh = fem::double0;
  int kflb = fem::int0;
  float xb = fem::float0;
  arr_1d<13, float> xfb(dim1(-6, 6), fem::fill0);
  double dshr = fem::double0;
  double dshz = fem::double0;
  float xe = fem::float0;
  float q2b = fem::float0;
  float tevb = fem::float0;
  float alsdum = fem::float0;
  float tevbsv = fem::float0;
  float b0 = fem::float0;
  arr_1d<13, float> wtap(dim1(-6, 6), fem::fill0);
  arr_1d<13, float> wtsf(dim1(-6, 6), fem::fill0);
  float wtapq = fem::float0;
  float wtsum = fem::float0;
  float xfbo = fem::float0;
  float q2ref = fem::float0;
  float wtran = fem::float0;
  int kfla = fem::int0;
  float z = fem::float0;
  float wtz = fem::float0;
  float rsoft = fem::float0;
  float alprat = fem::float0;
  float the2t = fem::float0;
  arr_1d<13, float> xfn(dim1(-6, 6), fem::fill0);
  float xfbn = fem::float0;
  float xa = fem::float0;
  arr_1d<13, float> xfa(dim1(-6, 6), fem::fill0);
  float xfan = fem::float0;
  float wtsfa = fem::float0;
  arr_1d<3, double> dq2(fem::fill0);
  double dplcm = fem::double0;
  int jr = fem::int0;
  int i = fem::int0;
  int ipo = fem::int0;
  int j = fem::int0;
  arr_1d<2, int> is(fem::fill0);
  arr_1d<3, double> dpc(fem::fill0);
  arr_1d<4, double> dpd(fem::fill0);
  int ikin = fem::int0;
  double dmsma = fem::double0;
  int it = fem::int0;
  float q2tim = fem::float0;
  double dms = fem::double0;
  double dpt2 = fem::double0;
  arr_1d<4, double> dpb(fem::fill0);
  double dbez = fem::double0;
  float the = fem::float0;
  int id1 = fem::int0;
  int id2 = fem::int0;
  int ir = fem::int0;
  arr_1d<5, float> robo(fem::fill0);
  float robot = fem::float0;
  //C
  //C...Generates spacelike parton showers.
  //C
  //C...Calculate maximum virtuality and check that evolution possible.
  ipus1 = ipu1;
  ipus2 = ipu2;
  isub = mint(1);
  q2e = vint(52);
  if (iset(isub) == 1) {
    q2e = q2e / parp(67);
  }
  else if (iset(isub) == 3 || iset(isub) == 4) {
    q2e = fem::pow2(pmas(23, 1));
    if (isub == 8 || isub == 76 || isub == 77) {
      q2e = fem::pow2(pmas(24, 1));
    }
  }
  tmax = fem::log(parp(67) * parp(63) * q2e / fem::pow2(parp(61)));
  if (parp(67) * q2e < fem::max(fem::pow2(parp(62)), 2.f * fem::pow2(
      parp(61))) || tmax < 0.2f) {
    return;
  }
  //C
  //C...Common constants and initial values. Save normal Lambda value.
  xe0 = 2.f * parp(65) / vint(1);
  alams = paru(111);
  paru(111) = parp(61);
  ns = n;
  statement_100:
  n = ns;
  FEM_DO_SAFE(jt, 1, 2) {
    kfls(jt) = mint(14 + jt);
    kfls(jt + 2) = kfls(jt);
    xs(jt) = vint(40 + jt);
    zs(jt) = 1.f;
    q2s(jt) = parp(67) * q2e;
    tevs(jt) = tmax;
    alam(jt) = parp(61);
    the2(jt) = 100.f;
    FEM_DO_SAFE(kfl, -6, 6) {
      xfs(jt, kfl) = xsfx(jt, kfl);
    }
  }
  dsh = fem::dble(vint(44));
  if (iset(isub) == 3 || iset(isub) == 4) {
    dsh = fem::dble(vint(26) * vint(2));
  }
  //C
  //C...Pick up leg with highest virtuality.
  statement_120:
  n++;
  jt = 1;
  if (n > ns + 1 && q2s(2) > q2s(1)) {
    jt = 2;
  }
  kflb = kfls(jt);
  xb = xs(jt);
  FEM_DO_SAFE(kfl, -6, 6) {
    xfb(kfl) = xfs(jt, kfl);
  }
  dshr = 2e0 * fem::sqrt(dsh);
  dshz = dsh / fem::dble(zs(jt));
  xe = fem::max(xe0, xb * (1.f / (1.f - parp(66)) - 1.f));
  if (xb + xe >= 0.999f) {
    q2b = 0.f;
    goto statement_220;
  }
  //C
  //C...Maximum Q2 without or with Q2 ordering. Effective Lambda and n_f.
  if (mstp(62) <= 1) {
    q2b = 0.5f * (1.f / zs(jt) + 1.f) * q2s(jt) + 0.5f * (1.f / zs(
      jt) - 1.f) * (q2s(3 - jt) - fem::sngl(dsh) + fem::sqrt(fem::pow2((
      fem::sngl(dsh) + q2s(1) + q2s(2))) + 8.f * q2s(1) * q2s(2) * zs(
      jt) / (1.f - zs(jt))));
    tevb = fem::log(parp(63) * q2b / fem::pow2(alam(jt)));
  }
  else {
    q2b = q2s(jt);
    tevb = tevs(jt);
  }
  alsdum = ulalps(cmn, parp(63) * q2b);
  tevb += 2.f * fem::log(alam(jt) / paru(117));
  tevbsv = tevb;
  alam(jt) = paru(117);
  b0 = (33.f - 2.f * mstu(118)) / 6.f;
  //C
  //C...Calculate Altarelli-Parisi and structure function weights.
  FEM_DO_SAFE(kfl, -6, 6) {
    wtap(kfl) = 0.f;
    wtsf(kfl) = 0.f;
  }
  if (kflb == 21) {
    wtapq = 16.f * (1.f - fem::sqrt(xb + xe)) / (3.f * fem::sqrt(xb));
    FEM_DO_SAFE(kfl, -mstp(54), mstp(54)) {
      if (kfl == 0) {
        wtap(kfl) = 6.f * fem::log((1.f - xb) / xe);
      }
      if (kfl != 0) {
        wtap(kfl) = wtapq;
      }
    }
  }
  else {
    wtap(0) = 0.5f * xb * (1.f / (xb + xe) - 1.f);
    wtap(kflb) = 8.f * fem::log((1.f - xb) * (xb + xe) / xe) / 3.f;
  }
  statement_160:
  wtsum = 0.f;
  if (kflb != 21) {
    xfbo = xfb(kflb);
  }
  if (kflb == 21) {
    xfbo = xfb(0);
  }
  //C***************************************************************
  //C**********ERROR HAS OCCURED HERE
  if (xfbo == 0.0f) {
    write(mstu(11), "(5x,'structure function has a zero point here')");
    write(mstu(11), "(5x,'xf(x,i=',i5,')=',f10.5)"), kflb, xfb(kflb);
    xfbo = 0.00001f;
  }
  //C****************************************************************
  FEM_DO_SAFE(kfl, -mstp(54), mstp(54)) {
    wtsf(kfl) = xfb(kfl) / xfbo;
    wtsum += wtap(kfl) * wtsf(kfl);
  }
  wtsum = fem::max(0.0001f, wtsum);
  //C
  //C...Choose new t: fix alpha_s, alpha_s(Q2), alpha_s(k_T2).
  statement_180:
  if (mstp(64) <= 0) {
    tevb += fem::log(rlu(cmn, 0)) * paru(2) / (paru(111) * wtsum);
  }
  else if (mstp(64) == 1) {
    tevb = tevb * fem::exp(fem::max(-100.f, fem::log(rlu(cmn, 0)) *
      b0 / wtsum));
  }
  else {
    tevb = tevb * fem::exp(fem::max(-100.f, fem::log(rlu(cmn, 0)) *
      b0 / (5.f * wtsum)));
  }
  statement_190:
  q2ref = fem::pow2(alam(jt)) * fem::exp(tevb);
  q2b = q2ref / parp(63);
  //C
  //C...Evolution ended or select flavour for branching parton.
  if (q2b < fem::pow2(parp(62))) {
    q2b = 0.f;
  }
  else {
    wtran = rlu(cmn, 0) * wtsum;
    kfla = -mstp(54) - 1;
    statement_200:
    kfla++;
    wtran = wtran - wtap(kfla) * wtsf(kfla);
    if (kfla < mstp(54) && wtran > 0.f) {
      goto statement_200;
    }
    if (kfla == 0) {
      kfla = 21;
    }
    //C
    //C...Choose z value and corrective weight.
    if (kflb == 21 && kfla == 21) {
      z = 1.f / (1.f + ((1.f - xb) / xb) * fem::pow((xe / (1.f - xb)),
        rlu(cmn, 0)));
      wtz = fem::pow2((1.f - z * (1.f - z)));
    }
    else if (kflb == 21) {
      z = xb / fem::pow2((1.f - rlu(cmn, 0) * (1.f - fem::sqrt(xb + xe))));
      wtz = 0.5f * (1.f + fem::pow2((1.f - z))) * fem::sqrt(z);
    }
    else if (kfla == 21) {
      z = xb * (1.f + rlu(cmn, 0) * (1.f / (xb + xe) - 1.f));
      wtz = 1.f - 2.f * z * (1.f - z);
    }
    else {
      z = 1.f - (1.f - xb) * fem::pow((xe / ((xb + xe) * (1.f - xb))),
        rlu(cmn, 0));
      wtz = 0.5f * (1.f + fem::pow2(z));
    }
    //C
    //C...Option with resummation of soft gluon emission as effective z shift.
    if (mstp(65) >= 1) {
      rsoft = 6.f;
      if (kflb != 21) {
        rsoft = 8.f / 3.f;
      }
      z = z * fem::pow((tevb / tevs(jt)), (rsoft * xe / ((xb + xe) * b0)));
      if (z <= xb) {
        goto statement_180;
      }
    }
    //C
    //C...Option with alpha_s(k_T2)Q2): demand k_T2 > cutoff, reweight.
    if (mstp(64) >= 2) {
      if ((1.f - z) * q2b < fem::pow2(parp(62))) {
        goto statement_180;
      }
      alprat = tevb / (tevb + fem::log(1.f - z));
      if (alprat < 5.f * rlu(cmn, 0)) {
        goto statement_180;
      }
      if (alprat > 5.f) {
        wtz = wtz * alprat / 5.f;
      }
    }
    //C
    //C...Option with angular ordering requirement.
    if (mstp(62) >= 3) {
      the2t = (4.f * fem::pow2(z) * q2b) / (vint(2) * (1.f - z) *
        fem::pow2(xb));
      if (the2t > the2(jt)) {
        goto statement_180;
      }
    }
    //C
    //C...Weighting with new structure functions.
    pystfu(cmn, mint(10 + jt), xb, q2ref, xfn, jt);
    if (kflb != 21) {
      xfbn = xfn(kflb);
    }
    if (kflb == 21) {
      xfbn = xfn(0);
    }
    if (xfbn < 1e-20f) {
      if (kfla == kflb) {
        tevb = tevbsv;
        wtap(kflb) = 0.f;
        goto statement_160;
      }
      else if (tevbsv - tevb > 0.2f) {
        tevb = 0.5f * (tevbsv + tevb);
        goto statement_190;
      }
      else {
        xfbn = 1e-10f;
      }
    }
    FEM_DO_SAFE(kfl, -mstp(54), mstp(54)) {
      xfb(kfl) = xfn(kfl);
    }
    xa = xb / z;
    pystfu(cmn, mint(10 + jt), xa, q2ref, xfa, jt);
    if (kfla != 21) {
      xfan = xfa(kfla);
    }
    if (kfla == 21) {
      xfan = xfa(0);
    }
    if (xfan < 1e-20f) {
      goto statement_160;
    }
    if (kfla != 21) {
      wtsfa = wtsf(kfla);
    }
    if (kfla == 21) {
      wtsfa = wtsf(0);
    }
    if (wtz * xfan / xfbn < rlu(cmn, 0) * wtsfa) {
      goto statement_160;
    }
  }
  //C
  //C...Define two hard scatterers in their CM-frame.
  statement_220:
  if (n == ns + 2) {
    dq2(jt) = fem::dble(q2b);
    dplcm = fem::dsqrt(fem::pow2((dsh + dq2(1) + dq2(2))) - 4e0 * dq2(
      1) * dq2(2)) / dshr;
    FEM_DO_SAFE(jr, 1, 2) {
      i = ns + jr;
      if (jr == 1) {
        ipo = ipus1;
      }
      if (jr == 2) {
        ipo = ipus2;
      }
      FEM_DO_SAFE(j, 1, 5) {
        k(i, j) = 0;
        p(i, j) = 0.f;
        v(i, j) = 0.f;
      }
      k(i, 1) = 14;
      k(i, 2) = kfls(jr + 2);
      k(i, 4) = ipo;
      k(i, 5) = ipo;
      p(i, 3) = fem::sngl(dplcm) * fem::pow((-1), (jr + 1));
      p(i, 4) = fem::sngl((dsh + dq2(3 - jr) - dq2(jr)) / dshr);
      p(i, 5) = -fem::sqrt(fem::sngl(dq2(jr)));
      k(ipo, 1) = 14;
      k(ipo, 3) = i;
      k(ipo, 4) = fem::mod(k(ipo, 4), mstu(5)) + mstu(5) * i;
      k(ipo, 5) = fem::mod(k(ipo, 5), mstu(5)) + mstu(5) * i;
    }
    //C
    //C...Find maximum allowed mass of timelike parton.
  }
  else if (n > ns + 2) {
    jr = 3 - jt;
    dq2(3) = fem::dble(q2b);
    dpc(1) = fem::dble(p(is(1), 4));
    dpc(2) = fem::dble(p(is(2), 4));
    dpc(3) = fem::dble(0.5f * (fem::abs(p(is(1), 3)) + fem::abs(p(is(2), 3))));
    dpd(1) = dsh + dq2(jr) + dq2(jt);
    dpd(2) = dshz + dq2(jr) + dq2(3);
    dpd(3) = fem::sqrt(fem::pow2(dpd(1)) - 4e0 * dq2(jr) * dq2(jt));
    dpd(4) = fem::sqrt(fem::pow2(dpd(2)) - 4e0 * dq2(jr) * dq2(3));
    ikin = 0;
    if (q2s(jr) >= fem::pow2((0.5f * parp(62))) && dpd(1) - dpd(
        3) >= 1e-10 * dpd(1)) {
      ikin = 1;
    }
    if (ikin == 0) {
      dmsma = (dq2(jt) / fem::dble(zs(jt)) - dq2(3)) * (dsh / (dsh +
        dq2(jt)) - dsh / (dshz + dq2(3)));
    }
    if (ikin == 1) {
      dmsma = (dpd(1) * dpd(2) - dpd(3) * dpd(4)) / (2.e0 * dq2(
        jr)) - dq2(jt) - dq2(3);
    }
    //C
    //C...Generate timelike parton shower (if required).
    it = n;
    FEM_DO_SAFE(j, 1, 5) {
      k(it, j) = 0;
      p(it, j) = 0.f;
      v(it, j) = 0.f;
    }
    k(it, 1) = 3;
    k(it, 2) = 21;
    if (kflb == 21 && kfls(jt + 2) != 21) {
      k(it, 2) = -kfls(jt + 2);
    }
    if (kflb != 21 && kfls(jt + 2) == 21) {
      k(it, 2) = kflb;
    }
    p(it, 5) = ulmass(cmn, k(it, 2));
    if (fem::sngl(dmsma) <= fem::pow2(p(it, 5))) {
      goto statement_100;
    }
    if (mstp(63) >= 1) {
      p(it, 4) = fem::sngl((dshz - dsh - fem::pow2(fem::dble(p(it,
        5)))) / dshr);
      p(it, 3) = fem::sqrt(fem::pow2(p(it, 4)) - fem::pow2(p(it, 5)));
      if (mstp(63) == 1) {
        q2tim = fem::sngl(dmsma);
      }
      else if (mstp(63) == 2) {
        q2tim = fem::min(fem::sngl(dmsma), parp(71) * q2s(jt));
      }
      else {
        //C'''Here remains to introduce angular ordering in first branching.
        q2tim = fem::sngl(dmsma);
      }
      lushow(cmn, it, 0, fem::sqrt(q2tim));
      if (n >= it + 1) {
        p(it, 5) = p(it + 1, 5);
      }
    }
    //C
    //C...Reconstruct kinematics of branching: timelike parton shower.
    dms = fem::dble(fem::pow2(p(it, 5)));
    if (ikin == 0) {
      dpt2 = (dmsma - dms) * (dshz + dq2(3)) / (dsh + dq2(jt));
    }
    if (ikin == 1) {
      dpt2 = (dmsma - dms) * (0.5e0 * dpd(1) * dpd(2) + 0.5e0 * dpd(
        3) * dpd(4) - dq2(jr) * (dq2(jt) + dq2(3) + dms)) / (4.e0 *
        dsh * fem::pow2(dpc(3)));
    }
    if (dpt2 < 0.e0) {
      goto statement_100;
    }
    dpb(1) = (0.5e0 * dpd(2) - dpc(jr) * (dshz + dq2(jr) - dq2(jt) -
      dms) / dshr) / dpc(3) - dpc(3);
    p(it, 1) = fem::sqrt(fem::sngl(dpt2));
    p(it, 3) = fem::sngl(dpb(1)) * fem::pow((-1), (jt + 1));
    p(it, 4) = fem::sngl((dshz - dsh - dms) / dshr);
    if (n >= it + 1) {
      dpb(1) = fem::sqrt(fem::pow2(dpb(1)) + dpt2);
      dpb(2) = fem::sqrt(fem::pow2(dpb(1)) + dms);
      dpb(3) = fem::dble(p(it + 1, 3));
      dpb(4) = fem::sqrt(fem::pow2(dpb(3)) + dms);
      dbez = (dpb(4) * dpb(1) - dpb(3) * dpb(2)) / (dpb(4) * dpb(2) -
        dpb(3) * dpb(1));
      ludbrb(it + 1, n, 0.f, 0.f, 0e0, 0e0, dbez);
      the = ulangl(cmn, p(it, 3), p(it, 1));
      ludbrb(it + 1, n, the, 0.f, 0e0, 0e0, 0e0);
    }
    //C
    //C...Reconstruct kinematics of branching: spacelike parton.
    FEM_DO_SAFE(j, 1, 5) {
      k(n + 1, j) = 0;
      p(n + 1, j) = 0.f;
      v(n + 1, j) = 0.f;
    }
    k(n + 1, 1) = 14;
    k(n + 1, 2) = kflb;
    p(n + 1, 1) = p(it, 1);
    p(n + 1, 3) = p(it, 3) + p(is(jt), 3);
    p(n + 1, 4) = p(it, 4) + p(is(jt), 4);
    p(n + 1, 5) = -fem::sqrt(fem::sngl(dq2(3)));
    //C
    //C...Define colour flow of branching.
    k(is(jt), 3) = n + 1;
    k(it, 3) = n + 1;
    id1 = it;
    if ((k(n + 1, 2) > 0 && k(n + 1, 2) != 21 && k(id1, 2) > 0 && k(id1,
        2) != 21) || (k(n + 1, 2) < 0 && k(id1, 2) == 21) || (k(n + 1,
        2) == 21 && k(id1, 2) == 21 && rlu(cmn, 0) > 0.5f) || (k(n + 1,
        2) == 21 && k(id1, 2) < 0)) {
      id1 = is(jt);
    }
    id2 = it + is(jt) - id1;
    k(n + 1, 4) += id1;
    k(n + 1, 5) += id2;
    k(id1, 4) += mstu(5) * (n + 1);
    k(id1, 5) += mstu(5) * id2;
    k(id2, 4) += mstu(5) * id1;
    k(id2, 5) += mstu(5) * (n + 1);
    n++;
    //C
    //C...Boost to new CM-frame.
    ludbrb(ns + 1, n, 0.f, 0.f, -fem::dble((p(n, 1) + p(is(jr), 1)) / (p(n,
      4) + p(is(jr), 4))), 0e0, -fem::dble((p(n, 3) + p(is(jr), 3)) / (p(n,
      4) + p(is(jr), 4))));
    ir = n + (jt - 1) * (is(1) - n);
    ludbrb(ns + 1, n, -ulangl(cmn, p(ir, 3), p(ir, 1)), paru(2) * rlu(cmn,
      0), 0e0, 0e0, 0e0);
  }
  //C
  //C...Save quantities, loop back.
  is(jt) = n;
  q2s(jt) = q2b;
  dq2(jt) = fem::dble(q2b);
  if (mstp(62) >= 3) {
    the2(jt) = the2t;
  }
  dsh = dshz;
  if (q2b >= fem::pow2((0.5f * parp(62)))) {
    kfls(jt + 2) = kfls(jt);
    kfls(jt) = kfla;
    xs(jt) = xa;
    zs(jt) = z;
    FEM_DO_SAFE(kfl, -6, 6) {
      xfs(jt, kfl) = xfa(kfl);
    }
    tevs(jt) = tevb;
  }
  else {
    if (jt == 1) {
      ipu1 = n;
    }
    if (jt == 2) {
      ipu2 = n;
    }
  }
  if (n > mstu(4) - mstu(32) - 10) {
    luerrm(cmn, 11, "(PYSSPA:) no more memory left in LUJETS");
    if (mstu(21) >= 1) {
      n = ns;
    }
    if (mstu(21) >= 1) {
      return;
    }
  }
  if (fem::max(q2s(1), q2s(2)) >= fem::pow2((0.5f * parp(62))) || n <= ns + 1) {
    goto statement_120;
  }
  //C
  //C...Boost hard scattering partons to frame of shower initiators.
  FEM_DO_SAFE(j, 1, 3) {
    robo(j + 2) = (p(ns + 1, j) + p(ns + 2, j)) / (p(ns + 1, 4) + p(ns + 2, 4));
  }
  FEM_DO_SAFE(j, 1, 5) {
    p(n + 2, j) = p(ns + 1, j);
  }
  robot = fem::pow2(robo(3)) + fem::pow2(robo(4)) + fem::pow2(robo(5));
  if (robot >= 0.999999f) {
    robot = 1.00001f * fem::sqrt(robot);
    robo(3) = robo(3) / robot;
    robo(4) = robo(4) / robot;
    robo(5) = robo(5) / robot;
  }
  ludbrb(n + 2, n + 2, 0.f, 0.f, -fem::dble(robo(3)), -fem::dble(robo(4)),
    -fem::dble(robo(5)));
  robo(2) = ulangl(cmn, p(n + 2, 1), p(n + 2, 2));
  robo(1) = ulangl(cmn, p(n + 2, 3), fem::sqrt(fem::pow2(p(n + 2,
    1)) + fem::pow2(p(n + 2, 2))));
  ludbrb(mint(83) + 5, ns, robo(1), robo(2), fem::dble(robo(3)),
    fem::dble(robo(4)), fem::dble(robo(5)));
  //C
  //C...Store user information. Reset Lambda value.
  k(ipu1, 3) = mint(83) + 3;
  k(ipu2, 3) = mint(83) + 4;
  FEM_DO_SAFE(jt, 1, 2) {
    mint(12 + jt) = kfls(jt);
    vint(140 + jt) = xs(jt);
  }
  paru(111) = alams;
  //C
}

//C
//C*********************************************************************
//C
void
pyremn(
  common& cmn,
  int const& ipu1,
  int const& ipu2)
{
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_cref<float, 2> pphi(cmn.pphi, dimension(300, 15));
  arr_cref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_cref<float, 2> pthi(cmn.pthi, dimension(300, 15));
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  //
  int jt = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  int isub = fem::int0;
  int ilep = fem::int0;
  int iq = fem::int0;
  int ip = fem::int0;
  int ilepr = fem::int0;
  int ns = fem::int0;
  int ipu = fem::int0;
  float shs = fem::float0;
  arr_1d<6, float> pms(fem::fill0);
  float rpt1 = fem::float0;
  float rpt2 = fem::float0;
  float ssw2 = fem::float0;
  float pt = fem::float0;
  float phi = fem::float0;
  int jpt = fem::int0;
  float ptgs = fem::float0;
  int iint = fem::int0;
  float pkcsq = fem::float0;
  float q2 = fem::float0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  float shr = fem::float0;
  arr_1d<5, float> robo(fem::fill0);
  int nmax = fem::int0;
  float peh = fem::float0;
  float pzh = fem::float0;
  float shh = fem::float0;
  float pmmin = fem::float0;
  float pei = fem::float0;
  float pzi = fem::float0;
  arr_1d<2, int> kflch(fem::fill0);
  arr_1d<2, int> kflsp(fem::fill0);
  arr_1d<2, int> is(fem::fill0);
  int kfls = fem::int0;
  int imb = fem::int0;
  float chik = fem::float0;
  arr_1d<2, float> chi(fem::fill0);
  float cut = fem::float0;
  float cutr = fem::float0;
  float chir = fem::float0;
  float pe = fem::float0;
  float pz = fem::float0;
  float pw1 = fem::float0;
  float pef = fem::float0;
  float pzf = fem::float0;
  float pt2 = fem::float0;
  float phipt = fem::float0;
  float rqp = fem::float0;
  float sinth = fem::float0;
  float betax = fem::float0;
  float pem = fem::float0;
  float pzm = fem::float0;
  float betaz = fem::float0;
  //C
  //C...Adds on target remnants (one or two from each side) and
  //C...includes primordial kT.
  //C...COMMON BLOCK FROM HIJING
  //C
  //C...Special case for lepton-lepton interaction.
  if (mint(43) == 1) {
    FEM_DO_SAFE(jt, 1, 2) {
      i = mint(83) + jt + 2;
      k(i, 1) = 21;
      k(i, 2) = k(i - 2, 2);
      k(i, 3) = i - 2;
      FEM_DO_SAFE(j, 1, 5) {
        p(i, j) = p(i - 2, j);
      }
    }
  }
  //C
  //C...Find event type, set pointers.
  if (ipu1 == 0 && ipu2 == 0) {
    return;
  }
  isub = mint(1);
  ilep = 0;
  if (ipu1 == 0) {
    ilep = 1;
  }
  if (ipu2 == 0) {
    ilep = 2;
  }
  if (isub == 95) {
    ilep = -1;
  }
  if (ilep == 1) {
    iq = mint(84) + 1;
  }
  if (ilep == 2) {
    iq = mint(84) + 2;
  }
  ip = fem::max(ipu1, ipu2);
  ilepr = mint(83) + 5 - ilep;
  ns = n;
  //C
  //C...Define initial partons, including primordial kT.
  statement_110:
  FEM_DO_SAFE(jt, 1, 2) {
    i = mint(83) + jt + 2;
    if (jt == 1) {
      ipu = ipu1;
    }
    if (jt == 2) {
      ipu = ipu2;
    }
    k(i, 1) = 21;
    k(i, 3) = i - 2;
    if (isub == 95) {
      k(i, 2) = 21;
      shs = 0.f;
    }
    else if (mint(40 + jt) == 1 && ipu != 0) {
      k(i, 2) = k(ipu, 2);
      p(i, 5) = p(ipu, 5);
      p(i, 1) = 0.f;
      p(i, 2) = 0.f;
      pms(jt) = fem::pow2(p(i, 5));
    }
    else if (ipu != 0) {
      k(i, 2) = k(ipu, 2);
      p(i, 5) = p(ipu, 5);
      //C...No primordial kT or chosen according to truncated Gaussian or
      //C...exponential.
      //C
      //C     X.N. Wang (7.22.97)
      //C
      rpt1 = 0.0f;
      rpt2 = 0.0f;
      ssw2 = fem::pow2((pphi(ihnt2(11), 4) + pthi(ihnt2(12), 4))) -
        fem::pow2((pphi(ihnt2(11), 1) + pthi(ihnt2(12), 1))) -
        fem::pow2((pphi(ihnt2(11), 2) + pthi(ihnt2(12), 2))) -
        fem::pow2((pphi(ihnt2(11), 3) + pthi(ihnt2(12), 3)));
      //C
      //C********this is s of the current NN collision
      if (ssw2 <= 4.0f * fem::pow2(parp(93))) {
        goto statement_1211;
      }
      //C
      if (ihpr2(5) <= 0) {
        statement_120:
        if (mstp(91) <= 0) {
          pt = 0.f;
        }
        else if (mstp(91) == 1) {
          pt = parp(91) * fem::sqrt(-fem::log(rlu(cmn, 0)));
        }
        else {
          rpt1 = rlu(cmn, 0);
          rpt2 = rlu(cmn, 0);
          pt = -parp(92) * fem::log(rpt1 * rpt2);
        }
        if (pt > parp(93)) {
          goto statement_120;
        }
        phi = paru(2) * rlu(cmn, 0);
        rpt1 = pt * fem::cos(phi);
        rpt2 = pt * fem::sin(phi);
      }
      else if (ihpr2(5) == 1) {
        if (jt == 1) {
          jpt = nfp(ihnt2(11), 11);
        }
        if (jt == 2) {
          jpt = nft(ihnt2(12), 11);
        }
        statement_1205:
        ptgs = parp(91) * fem::sqrt(-fem::log(rlu(cmn, 0)));
        if (ptgs > parp(93)) {
          goto statement_1205;
        }
        phi = 2.0f * hipr1(40) * rlu(cmn, 0);
        rpt1 = ptgs * fem::cos(phi);
        rpt2 = ptgs * fem::sin(phi);
        FEM_DO_SAFE(iint, 1, jpt - 1) {
          pkcsq = parp(91) * fem::sqrt(-fem::log(rlu(cmn, 0)));
          phi = 2.0f * hipr1(40) * rlu(cmn, 0);
          rpt1 += pkcsq * fem::cos(phi);
          rpt2 += pkcsq * fem::sin(phi);
        }
        if (fem::pow2(rpt1) + fem::pow2(rpt2) >= ssw2 / 4.0f) {
          goto statement_1205;
        }
      }
      //C     X.N. Wang
      //C                     ********When initial interaction among soft partons is
      //C                             assumed the primordial pt comes from the sum of
      //C                             pt of JPT-1 number of initial interaction, JPT
      //C                             is the number of interaction including present
      //C                             one that nucleon hassuffered
      statement_1211:
      p(i, 1) = rpt1;
      p(i, 2) = rpt2;
      pms(jt) = fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) + fem::pow2(p(i, 2));
    }
    else {
      k(i, 2) = k(iq, 2);
      q2 = vint(52);
      p(i, 5) = -fem::sqrt(q2);
      pms(jt) = -q2;
      shs = (1.f - vint(43 - jt)) * q2 / vint(43 - jt) + fem::pow2(
        vint(5 - jt));
    }
  }
  //C
  //C...Kinematics construction for initial partons.
  i1 = mint(83) + 3;
  i2 = mint(83) + 4;
  if (ilep == 0) {
    shs = vint(141) * vint(142) * vint(2) + fem::pow2((p(i1, 1) + p(i2,
      1))) + fem::pow2((p(i1, 2) + p(i2, 2)));
  }
  shr = fem::sqrt(fem::max(0.f, shs));
  if (ilep == 0) {
    if (fem::pow2((shs - pms(1) - pms(2))) - 4.f * pms(1) * pms(2) <= 0.f) {
      goto statement_110;
    }
    p(i1, 4) = 0.5f * (shr + (pms(1) - pms(2)) / shr);
    p(i1, 3) = fem::sqrt(fem::max(0.f, fem::pow2(p(i1, 4)) - pms(1)));
    p(i2, 4) = shr - p(i1, 4);
    p(i2, 3) = -p(i1, 3);
  }
  else if (ilep == 1) {
    p(i1, 4) = p(iq, 4);
    p(i1, 3) = p(iq, 3);
    p(i2, 4) = p(ip, 4);
    p(i2, 3) = p(ip, 3);
  }
  else if (ilep == 2) {
    p(i1, 4) = p(ip, 4);
    p(i1, 3) = p(ip, 3);
    p(i2, 4) = p(iq, 4);
    p(i2, 3) = p(iq, 3);
  }
  if (mint(43) == 1) {
    return;
  }
  //C
  //C...Transform partons to overall CM-frame (not for leptoproduction).
  if (ilep == 0) {
    robo(3) = (p(i1, 1) + p(i2, 1)) / shr;
    robo(4) = (p(i1, 2) + p(i2, 2)) / shr;
    ludbrb(i1, i2, 0.f, 0.f, -fem::dble(robo(3)), -fem::dble(robo(4)), 0e0);
    robo(2) = ulangl(cmn, p(i1, 1), p(i1, 2));
    ludbrb(i1, i2, 0.f, -robo(2), 0e0, 0e0, 0e0);
    robo(1) = ulangl(cmn, p(i1, 3), p(i1, 1));
    ludbrb(i1, i2, -robo(1), 0.f, 0e0, 0e0, 0e0);
    nmax = fem::max(mint(52), ipu1, ipu2);
    ludbrb(i1, nmax, robo(1), robo(2), fem::dble(robo(3)), fem::dble(robo(4)),
      0e0);
    robo(5) = fem::max(-0.999999f, fem::min(0.999999f, (vint(141) -
      vint(142)) / (vint(141) + vint(142))));
    ludbrb(i1, nmax, 0.f, 0.f, 0e0, 0e0, fem::dble(robo(5)));
  }
  //C
  //C...Check invariant mass of remnant system:
  //C...hadronic events or leptoproduction.
  if (ilep <= 0) {
    if (mstp(81) <= 0 || mstp(82) <= 0 || isub == 95) {
      vint(151) = 0.f;
      vint(152) = 0.f;
    }
    peh = p(i1, 4) + p(i2, 4) + 0.5f * vint(1) * (vint(151) + vint(152));
    pzh = p(i1, 3) + p(i2, 3) + 0.5f * vint(1) * (vint(151) - vint(152));
    shh = fem::pow2((vint(1) - peh)) - fem::pow2((p(i1, 1) + p(i2,
      1))) - fem::pow2((p(i1, 2) + p(i2, 2))) - fem::pow2(pzh);
    pmmin = p(mint(83) + 1, 5) + p(mint(83) + 2, 5) + ulmass(cmn, k(i1,
      2)) + ulmass(cmn, k(i2, 2));
    if (shr >= vint(1) || shh <= fem::pow2((pmmin + parp(111)))) {
      mint(51) = 1;
      return;
    }
    shr = fem::sqrt(shh + fem::pow2((p(i1, 1) + p(i2, 1))) + fem::pow2((p(i1,
      2) + p(i2, 2))));
  }
  else {
    pei = p(iq, 4) + p(ip, 4);
    pzi = p(iq, 3) + p(ip, 3);
    pms(ilep) = fem::max(0.f, fem::pow2(pei) - fem::pow2(pzi));
    pmmin = p(ilepr - 2, 5) + ulmass(cmn, k(ilepr, 2)) + fem::sqrt(pms(ilep));
    if (shr <= pmmin + parp(111)) {
      mint(51) = 1;
      return;
    }
  }
  //C
  //C...Subdivide remnant if necessary, store first parton.
  statement_140:
  i = ns;
  FEM_DO_SAFE(jt, 1, 2) {
    if (jt == ilep) {
      goto statement_190;
    }
    if (jt == 1) {
      ipu = ipu1;
    }
    if (jt == 2) {
      ipu = ipu2;
    }
    pyspli(cmn, mint(10 + jt), mint(12 + jt), kflch(jt), kflsp(jt));
    i++;
    is(jt) = i;
    FEM_DO_SAFE(j, 1, 5) {
      k(i, j) = 0;
      p(i, j) = 0.f;
      v(i, j) = 0.f;
    }
    k(i, 1) = 3;
    k(i, 2) = kflsp(jt);
    k(i, 3) = mint(83) + jt;
    p(i, 5) = ulmass(cmn, k(i, 2));
    //C
    //C...First parton colour connections and transverse mass.
    kfls = (3 - kchg(lucomp(cmn, kflsp(jt)), 2) * fem::isign(1, kflsp(jt))) / 2;
    k(i, kfls + 3) = ipu;
    k(ipu, 6 - kfls) = fem::mod(k(ipu, 6 - kfls), mstu(5)) + mstu(5) * i;
    if (kflch(jt) == 0) {
      p(i, 1) = -p(mint(83) + jt + 2, 1);
      p(i, 2) = -p(mint(83) + jt + 2, 2);
      pms(jt) = fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) + fem::pow2(p(i, 2));
      //C
      //C...When extra remnant parton or hadron: find relative pT, store.
    }
    else {
      luptdi(cmn, 1, p(i, 1), p(i, 2));
      pms(jt + 2) = fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) + fem::pow2(p(i,
        2));
      i++;
      FEM_DO_SAFE(j, 1, 5) {
        k(i, j) = 0;
        p(i, j) = 0.f;
        v(i, j) = 0.f;
      }
      k(i, 1) = 1;
      k(i, 2) = kflch(jt);
      k(i, 3) = mint(83) + jt;
      p(i, 5) = ulmass(cmn, k(i, 2));
      p(i, 1) = -p(mint(83) + jt + 2, 1) - p(i - 1, 1);
      p(i, 2) = -p(mint(83) + jt + 2, 2) - p(i - 1, 2);
      pms(jt + 4) = fem::pow2(p(i, 5)) + fem::pow2(p(i, 1)) + fem::pow2(p(i,
        2));
      //C...Relative distribution of energy for particle into two jets.
      imb = 1;
      if (fem::mod(mint(10 + jt) / 1000, 10) != 0) {
        imb = 2;
      }
      if (fem::iabs(kflch(jt)) <= 10 || kflch(jt) == 21) {
        chik = parp(92 + 2 * imb);
        if (mstp(92) <= 1) {
          if (imb == 1) {
            chi(jt) = rlu(cmn, 0);
          }
          if (imb == 2) {
            chi(jt) = 1.f - fem::sqrt(rlu(cmn, 0));
          }
        }
        else if (mstp(92) == 2) {
          chi(jt) = 1.f - fem::pow(rlu(cmn, 0), (1.f / (1.f + chik)));
        }
        else if (mstp(92) == 3) {
          cut = 2.f * 0.3f / vint(1);
          statement_170:
          chi(jt) = fem::pow2(rlu(cmn, 0));
          if (fem::pow((fem::pow2(chi(jt)) / (fem::pow2(chi(jt)) +
              fem::pow2(cut))), 0.25f) * fem::pow((1.f - chi(jt)),
              chik) < rlu(cmn, 0)) {
            goto statement_170;
          }
        }
        else {
          cut = 2.f * 0.3f / vint(1);
          cutr = (1.f + fem::sqrt(1.f + fem::pow2(cut))) / cut;
          statement_180:
          chir = cut * fem::pow(cutr, rlu(cmn, 0));
          chi(jt) = (fem::pow2(chir) - fem::pow2(cut)) / (2.f * chir);
          if (fem::pow((1.f - chi(jt)), chik) < rlu(cmn, 0)) {
            goto statement_180;
          }
        }
        //C...Relative distribution of energy for particle into jet plus particle.
      }
      else {
        if (mstp(92) <= 1) {
          if (imb == 1) {
            chi(jt) = rlu(cmn, 0);
          }
          if (imb == 2) {
            chi(jt) = 1.f - fem::sqrt(rlu(cmn, 0));
          }
        }
        else {
          chi(jt) = 1.f - fem::pow(rlu(cmn, 0), (1.f / (1.f + parp(
            93 + 2 * imb))));
        }
        if (fem::mod(kflch(jt) / 1000, 10) != 0) {
          chi(jt) = 1.f - chi(jt);
        }
      }
      pms(jt) = pms(jt + 4) / chi(jt) + pms(jt + 2) / (1.f - chi(jt));
      kfls = kchg(lucomp(cmn, kflch(jt)), 2) * fem::isign(1, kflch(jt));
      if (kfls != 0) {
        k(i, 1) = 3;
        kfls = (3 - kfls) / 2;
        k(i, kfls + 3) = ipu;
        k(ipu, 6 - kfls) = fem::mod(k(ipu, 6 - kfls), mstu(5)) + mstu(5) * i;
      }
    }
    statement_190:;
  }
  if (shr <= fem::sqrt(pms(1)) + fem::sqrt(pms(2))) {
    goto statement_140;
  }
  n = i;
  //C
  //C...Reconstruct kinematics of remnants.
  FEM_DO_SAFE(jt, 1, 2) {
    if (jt == ilep) {
      goto statement_200;
    }
    pe = 0.5f * (shr + (pms(jt) - pms(3 - jt)) / shr);
    pz = fem::sqrt(fem::pow2(pe) - pms(jt));
    if (kflch(jt) == 0) {
      p(is(jt), 4) = pe;
      p(is(jt), 3) = pz * fem::pow((-1), (jt - 1));
    }
    else {
      pw1 = chi(jt) * (pe + pz);
      p(is(jt) + 1, 4) = 0.5f * (pw1 + pms(jt + 4) / pw1);
      p(is(jt) + 1, 3) = 0.5f * (pw1 - pms(jt + 4) / pw1) * fem::pow((-1),
        (jt - 1));
      p(is(jt), 4) = pe - p(is(jt) + 1, 4);
      p(is(jt), 3) = pz * fem::pow((-1), (jt - 1)) - p(is(jt) + 1, 3);
    }
    statement_200:;
  }
  //C
  //C...Hadronic events: boost remnants to correct longitudinal frame.
  if (ilep <= 0) {
    ludbrb(ns + 1, n, 0.f, 0.f, 0e0, 0e0, -fem::dble(pzh / (vint(1) - peh)));
    //C...Leptoproduction events: boost colliding subsystem.
  }
  else {
    nmax = fem::max(ip, mint(52));
    pef = shr - pe;
    pzf = pz * fem::pow((-1), (ilep - 1));
    pt2 = fem::pow2(p(ilepr, 1)) + fem::pow2(p(ilepr, 2));
    phipt = ulangl(cmn, p(ilepr, 1), p(ilepr, 2));
    ludbrb(mint(84) + 1, nmax, 0.f, -phipt, 0e0, 0e0, 0e0);
    rqp = p(iq, 3) * (pt2 + fem::pow2(pei)) - p(iq, 4) * pei * pzi;
    sinth = p(iq, 4) * fem::sqrt(pt2 * (pt2 + fem::pow2(pei)) / (
      fem::pow2(rqp) + pt2 * fem::pow2(p(iq, 4)) * fem::pow2(pzi))) *
      fem::sign(1.f, -rqp);
    ludbrb(mint(84) + 1, nmax, fem::asin(sinth), 0.f, 0e0, 0e0, 0e0);
    betax = (-pei * pzi * sinth + fem::sqrt(pt2 * (pt2 + fem::pow2(
      pei) - fem::pow2((pzi * sinth))))) / (pt2 + fem::pow2(pei));
    ludbrb(mint(84) + 1, nmax, 0.f, 0.f, fem::dble(betax), 0e0, 0e0);
    ludbrb(mint(84) + 1, nmax, 0.f, phipt, 0e0, 0e0, 0e0);
    pem = p(iq, 4) + p(ip, 4);
    pzm = p(iq, 3) + p(ip, 3);
    betaz = (-pem * pzm + pzf * fem::sqrt(fem::pow2(pzf) + fem::pow2(
      pem) - fem::pow2(pzm))) / (fem::pow2(pzf) + fem::pow2(pem));
    ludbrb(mint(84) + 1, nmax, 0.f, 0.f, 0e0, 0e0, fem::dble(betaz));
    ludbrb(i1, i2, fem::asin(sinth), 0.f, fem::dble(betax), 0e0, 0e0);
    ludbrb(i1, i2, 0.f, phipt, 0e0, 0e0, fem::dble(betaz));
  }
  //C
}

//C
//C*********************************************************************
//C
void
pyresd(
  common& cmn)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_cref<float> paru(cmn.paru, dimension(200));
  arr_cref<int, 2> kchg(cmn.kchg, dimension(500, 3));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_cref<int, 2> kfdp(cmn.kfdp, dimension(2000, 5));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_cref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  //
  int i1 = fem::int0;
  int i2 = fem::int0;
  int i3 = fem::int0;
  int i4 = fem::int0;
  int i5 = fem::int0;
  int i6 = fem::int0;
  arr_2d<6, 6, std::complex<float> > ha(fem::fill0);
  arr_2d<6, 6, std::complex<float> > hc(fem::fill0);
  double dt = fem::double0;
  double du = fem::double0;
  double d34 = fem::double0;
  double d56 = fem::double0;
  int isub = fem::int0;
  float sh = fem::float0;
  arr_2d<10, 6, int> iref(fem::fill0);
  int np = fem::int0;
  int ip = fem::int0;
  int ninh = fem::int0;
  int jtmax = fem::int0;
  int jt = fem::int0;
  arr_1d<2, int> kdcy(fem::fill0);
  arr_1d<2, int> kfl1(fem::fill0);
  arr_1d<2, int> kfl2(fem::fill0);
  arr_1d<2, int> nsd(fem::fill0);
  int id = fem::int0;
  int kfa = fem::int0;
  arr_1d<41, float> wdtp(dim1(0, 40), fem::fill0);
  arr_2d<41, 6, float> wdte(dim1(0, 40).dim2(0, 5), fem::fill0);
  int ipm = fem::int0;
  int i12 = fem::int0;
  float rkfl = fem::float0;
  int i = fem::int0;
  int idc = fem::int0;
  float pid5 = fem::float0;
  arr_1d<2, float> cthe(fem::fill0);
  arr_1d<2, float> phi(fem::fill0);
  arr_1d<6, int> ilin(fem::fill0);
  int imin = fem::int0;
  int imax = fem::int0;
  int iord = fem::int0;
  float xw = fem::float0;
  int j = fem::int0;
  arr_2d<6, 4, float> coup(fem::fill0);
  float sqmz = fem::float0;
  float gzmz = fem::float0;
  float sqmw = fem::float0;
  float gzmw = fem::float0;
  float sqmzp = fem::float0;
  float gzmzp = fem::float0;
  float therr = fem::float0;
  float phirr = fem::float0;
  arr_2d<6, 4, float> pk(fem::fill0);
  arr_2d<6, 6, float> pkk(fem::fill0);
  float wt = fem::float0;
  float wtmax = fem::float0;
  float ei = fem::float0;
  float ai = fem::float0;
  float vi = fem::float0;
  float ef = fem::float0;
  float af = fem::float0;
  float vf = fem::float0;
  float gg = fem::float0;
  float gz = fem::float0;
  float zz = fem::float0;
  float asym = fem::float0;
  float s34 = fem::float0;
  float s56 = fem::float0;
  float ti = fem::float0;
  float ui = fem::float0;
  float cawz = fem::float0;
  float cbwz = fem::float0;
  float cdww = fem::float0;
  float caww = fem::float0;
  float cbww = fem::float0;
  float ccww = fem::float0;
  float api = fem::float0;
  float vpi = fem::float0;
  float apf = fem::float0;
  float vpf = fem::float0;
  float gzp = fem::float0;
  float zzp = fem::float0;
  float zpzp = fem::float0;
  int idoc = fem::int0;
  //C
  //C...Allows resonances to decay (including parton showers for hadronic
  //C...channels).
  //C
  //C...The F, Xi and Xj functions of Gunion and Kunszt
  //C...(Phys. Rev. D33, 665, plus errata from the authors).
  fgk(i1, i2, i3, i4, i5, i6) = 4.f * ha(i1, i3) * hc(i2, i6) * (ha(i1,
    i5) * hc(i1, i4) + ha(i3, i5) * hc(i3, i4));
  digk(dt, du) = -4.e0 * d34 * d56 + dt * (3.e0 * dt + 4.e0 * du) +
    fem::pow2(dt) * (dt * du / (d34 * d56) - 2.e0 * (1.e0 / d34 +
    1.e0 / d56) * (dt + du) + 2.e0 * (d34 / d56 + d56 / d34));
  djgk(dt, du) = 8.e0 * fem::pow2((d34 + d56)) - 8.e0 * (d34 + d56) * (
    dt + du) - 6.e0 * dt * du - 2.e0 * dt * du * (dt * du / (d34 *
    d56) - 2.e0 * (1.e0 / d34 + 1.e0 / d56) * (dt + du) + 2.e0 * (d34 /
    d56 + d56 / d34));
  //C
  //C...Define initial two objects, initialize loop.
  isub = mint(1);
  sh = vint(44);
  iref(1, 5) = 0;
  iref(1, 6) = 0;
  if (iset(isub) == 1 || iset(isub) == 3) {
    iref(1, 1) = mint(84) + 2 + iset(isub);
    iref(1, 2) = 0;
    iref(1, 3) = mint(83) + 6 + iset(isub);
    iref(1, 4) = 0;
  }
  else if (iset(isub) == 2 || iset(isub) == 4) {
    iref(1, 1) = mint(84) + 1 + iset(isub);
    iref(1, 2) = mint(84) + 2 + iset(isub);
    iref(1, 3) = mint(83) + 5 + iset(isub);
    iref(1, 4) = mint(83) + 6 + iset(isub);
  }
  np = 1;
  ip = 0;
  statement_100:
  ip++;
  ninh = 0;
  //C
  //C...Loop over one/two resonances; reset decay rates.
  jtmax = 2;
  if (ip == 1 && (iset(isub) == 1 || iset(isub) == 3)) {
    jtmax = 1;
  }
  FEM_DO_SAFE(jt, 1, jtmax) {
    kdcy(jt) = 0;
    kfl1(jt) = 0;
    kfl2(jt) = 0;
    nsd(jt) = iref(ip, jt);
    id = iref(ip, jt);
    if (id == 0) {
      goto statement_140;
    }
    kfa = fem::iabs(k(id, 2));
    if (kfa < 23 || kfa > 40) {
      goto statement_140;
    }
    if (mdcy(kfa, 1) != 0) {
      if (isub == 1 || isub == 141) {
        mint(61) = 1;
      }
      pywidt(cmn, kfa, p(id, 5), wdtp, wdte);
      if (kchg(kfa, 3) == 0) {
        ipm = 2;
      }
      else {
        ipm = (5 + fem::isign(1, k(id, 2))) / 2;
      }
      if (jtmax == 1 || fem::iabs(k(iref(ip, 1), 2)) != fem::iabs(k(iref(ip,
          2), 2))) {
        i12 = 4;
      }
      else {
        if (jt == 1) {
          i12 = fem::fint(4.5f + rlu(cmn, 0));
        }
        i12 = 9 - i12;
      }
      rkfl = (wdte(0, 1) + wdte(0, ipm) + wdte(0, i12)) * rlu(cmn, 0);
      FEM_DO_SAFE(i, 1, mdcy(kfa, 3)) {
        idc = i + mdcy(kfa, 2) - 1;
        kfl1(jt) = kfdp(idc, 1) * fem::isign(1, k(id, 2));
        kfl2(jt) = kfdp(idc, 2) * fem::isign(1, k(id, 2));
        rkfl = rkfl - (wdte(i, 1) + wdte(i, ipm) + wdte(i, i12));
        if (rkfl <= 0.f) {
          goto statement_130;
        }
      }
      statement_130:;
    }
    //C
    //C...Summarize result on decay channel chosen.
    if ((kfa == 23 || kfa == 24) && kfl1(jt) == 0) {
      ninh++;
    }
    if (kfl1(jt) == 0) {
      goto statement_140;
    }
    kdcy(jt) = 2;
    if (fem::iabs(kfl1(jt)) <= 10 || kfl1(jt) == 21) {
      kdcy(jt) = 1;
    }
    if ((fem::iabs(kfl1(jt)) >= 23 && fem::iabs(kfl1(jt)) <= 25) || (
        fem::iabs(kfl1(jt)) == 37)) {
      kdcy(jt) = 3;
    }
    nsd(jt) = n;
    //C
    //C...Fill decay products, prepared for parton showers for quarks.
    //Clin-8/19/02 avoid actual argument in common blocks of LU2ENT:
    pid5 = p(id, 5);
    if (kdcy(jt) == 1) {
      //C        CALL LU2ENT(-(N+1),KFL1(JT),KFL2(JT),P(ID,5))
      lu2ent(cmn, -(n + 1), kfl1(jt), kfl2(jt), pid5);
    }
    else {
      //C        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
      lu2ent(cmn, n + 1, kfl1(jt), kfl2(jt), pid5);
    }
    //C
    if (jtmax == 1) {
      cthe(jt) = vint(13) + (vint(33) - vint(13) + vint(34) - vint(
        14)) * rlu(cmn, 0);
      if (cthe(jt) > vint(33)) {
        cthe(jt) += vint(14) - vint(33);
      }
      phi(jt) = vint(24);
    }
    else {
      cthe(jt) = 2.f * rlu(cmn, 0) - 1.f;
      phi(jt) = paru(2) * rlu(cmn, 0);
    }
    statement_140:;
  }
  if (mint(3) == 1 && ip == 1) {
    mint(25) = kfl1(1);
    mint(26) = kfl2(1);
  }
  if (jtmax == 1 && kdcy(1) == 0) {
    goto statement_530;
  }
  if (jtmax == 2 && kdcy(1) == 0 && kdcy(2) == 0) {
    goto statement_530;
  }
  if (mstp(45) <= 0 || iref(ip, 2) == 0 || ninh >= 1) {
    goto statement_500;
  }
  if (k(iref(1, 1), 2) == 25 && ip == 1) {
    goto statement_500;
  }
  if (k(iref(1, 1), 2) == 25 && kdcy(1) * kdcy(2) == 0) {
    goto statement_500;
  }
  //C
  //C...Order incoming partons and outgoing resonances.
  ilin(1) = mint(84) + 1;
  if (k(mint(84) + 1, 2) > 0) {
    ilin(1) = mint(84) + 2;
  }
  if (k(ilin(1), 2) == 21) {
    ilin(1) = 2 * mint(84) + 3 - ilin(1);
  }
  ilin(2) = 2 * mint(84) + 3 - ilin(1);
  imin = 1;
  if (iref(ip, 5) == 25) {
    imin = 3;
  }
  imax = 2;
  iord = 1;
  if (k(iref(ip, 1), 2) == 23) {
    iord = 2;
  }
  if (k(iref(ip, 1), 2) == 24 && k(iref(ip, 2), 2) ==  - 24) {
    iord = 2;
  }
  if (fem::iabs(k(iref(ip, iord), 2)) == 25) {
    iord = 3 - iord;
  }
  if (kdcy(iord) == 0) {
    iord = 3 - iord;
  }
  //C
  //C...Order decay products of resonances.
  FEM_DOSTEP(jt, iord, 3 - iord, 3 - 2 * iord) {
    if (kdcy(jt) == 0) {
      ilin(imax + 1) = nsd(jt);
      imax++;
    }
    else if (k(nsd(jt) + 1, 2) > 0) {
      ilin(imax + 1) = n + 2 * jt - 1;
      ilin(imax + 2) = n + 2 * jt;
      imax += 2;
      k(n + 2 * jt - 1, 2) = k(nsd(jt) + 1, 2);
      k(n + 2 * jt, 2) = k(nsd(jt) + 2, 2);
    }
    else {
      ilin(imax + 1) = n + 2 * jt;
      ilin(imax + 2) = n + 2 * jt - 1;
      imax += 2;
      k(n + 2 * jt - 1, 2) = k(nsd(jt) + 1, 2);
      k(n + 2 * jt, 2) = k(nsd(jt) + 2, 2);
    }
  }
  //C
  //C...Find charge, isospin, left- and righthanded couplings.
  xw = paru(102);
  FEM_DO_SAFE(i, imin, imax) {
    FEM_DO_SAFE(j, 1, 4) {
      coup(i, j) = 0.f;
    }
    kfa = fem::iabs(k(ilin(i), 2));
    if (kfa > 20) {
      goto statement_410;
    }
    coup(i, 1) = luchge(cmn, kfa) / 3.f;
    coup(i, 2) = fem::pow((-1), fem::mod(kfa, 2));
    coup(i, 4) = -2.f * coup(i, 1) * xw;
    coup(i, 3) = coup(i, 2) + coup(i, 4);
    statement_410:;
  }
  sqmz = fem::pow2(pmas(23, 1));
  gzmz = pmas(23, 1) * pmas(23, 2);
  sqmw = fem::pow2(pmas(24, 1));
  gzmw = pmas(24, 1) * pmas(24, 2);
  sqmzp = fem::pow2(pmas(32, 1));
  gzmzp = pmas(32, 1) * pmas(32, 2);
  //C
  //C...Select random angles; construct massless four-vectors.
  statement_420:
  FEM_DO_SAFE(i, n + 1, n + 4) {
    k(i, 1) = 1;
    FEM_DO_SAFE(j, 1, 5) {
      p(i, j) = 0.f;
    }
  }
  FEM_DO_SAFE(jt, 1, jtmax) {
    if (kdcy(jt) == 0) {
      goto statement_440;
    }
    id = iref(ip, jt);
    p(n + 2 * jt - 1, 3) = 0.5f * p(id, 5);
    p(n + 2 * jt - 1, 4) = 0.5f * p(id, 5);
    p(n + 2 * jt, 3) = -0.5f * p(id, 5);
    p(n + 2 * jt, 4) = 0.5f * p(id, 5);
    cthe(jt) = 2.f * rlu(cmn, 0) - 1.f;
    phi(jt) = paru(2) * rlu(cmn, 0);
    ludbrb(n + 2 * jt - 1, n + 2 * jt, fem::acos(cthe(jt)), phi(jt),
      fem::dble(p(id, 1) / p(id, 4)), fem::dble(p(id, 2) / p(id, 4)),
      fem::dble(p(id, 3) / p(id, 4)));
    statement_440:;
  }
  //C
  //C...Store incoming and outgoing momenta, with random rotation to
  //C...avoid accidental zeroes in HA expressions.
  FEM_DO_SAFE(i, 1, imax) {
    k(n + 4 + i, 1) = 1;
    p(n + 4 + i, 4) = fem::sqrt(fem::pow2(p(ilin(i), 1)) + fem::pow2(p(ilin(i),
      2)) + fem::pow2(p(ilin(i), 3)) + fem::pow2(p(ilin(i), 5)));
    p(n + 4 + i, 5) = p(ilin(i), 5);
    FEM_DO_SAFE(j, 1, 3) {
      p(n + 4 + i, j) = p(ilin(i), j);
    }
  }
  therr = fem::acos(2.f * rlu(cmn, 0) - 1.f);
  phirr = paru(2) * rlu(cmn, 0);
  ludbrb(n + 5, n + 4 + imax, therr, phirr, 0e0, 0e0, 0e0);
  FEM_DO_SAFE(i, 1, imax) {
    FEM_DO_SAFE(j, 1, 4) {
      pk(i, j) = p(n + 4 + i, j);
    }
  }
  //C
  //C...Calculate internal products.
  if (isub == 22 || isub == 23 || isub == 25) {
    FEM_DO_SAFE(i1, imin, imax - 1) {
      FEM_DO_SAFE(i2, i1 + 1, imax) {
        ha(i1, i2) = fem::sqrt((pk(i1, 4) - pk(i1, 3)) * (pk(i2, 4) + pk(i2,
          3)) / (1e-20f + fem::pow2(pk(i1, 1)) + fem::pow2(pk(i1,
          2)))) * fem::cmplx(pk(i1, 1), pk(i1, 2)) - fem::sqrt((pk(i1,
          4) + pk(i1, 3)) * (pk(i2, 4) - pk(i2, 3)) / (1e-20f +
          fem::pow2(pk(i2, 1)) + fem::pow2(pk(i2, 2)))) * fem::cmplx(pk(i2,
          1), pk(i2, 2));
        hc(i1, i2) = fem::conjg(ha(i1, i2));
        if (i1 <= 2) {
          ha(i1, i2) = fem::cmplx(0.f, 1.f) * ha(i1, i2);
        }
        if (i1 <= 2) {
          hc(i1, i2) = fem::cmplx(0.f, 1.f) * hc(i1, i2);
        }
        ha(i2, i1) = -ha(i1, i2);
        hc(i2, i1) = -hc(i1, i2);
      }
    }
  }
  FEM_DO_SAFE(i, 1, 2) {
    FEM_DO_SAFE(j, 1, 4) {
      pk(i, j) = -pk(i, j);
    }
  }
  FEM_DO_SAFE(i1, imin, imax - 1) {
    FEM_DO_SAFE(i2, i1 + 1, imax) {
      pkk(i1, i2) = 2.f * (pk(i1, 4) * pk(i2, 4) - pk(i1, 1) * pk(i2,
        1) - pk(i1, 2) * pk(i2, 2) - pk(i1, 3) * pk(i2, 3));
      pkk(i2, i1) = pkk(i1, i2);
    }
  }
  //C
  if (iref(ip, 5) == 25) {
    //C...Angular weight for H0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons
    wt = 16.f * pkk(3, 5) * pkk(4, 6);
    if (ip == 1) {
      wtmax = fem::pow2(sh);
    }
    if (ip >= 2) {
      wtmax = fem::pow4(p(iref(ip, 6), 5));
    }
    //C
  }
  else if (isub == 1) {
    if (kfa != 37) {
      //C...Angular weight for gamma*/Z0 -> 2 quarks/leptons
      ei = kchg(fem::iabs(mint(15)), 1) / 3.f;
      ai = fem::sign(1.f, ei + 0.1f);
      vi = ai - 4.f * ei * xw;
      ef = kchg(kfa, 1) / 3.f;
      af = fem::sign(1.f, ef + 0.1f);
      vf = af - 4.f * ef * xw;
      gg = 1.f;
      gz = 1.f / (8.f * xw * (1.f - xw)) * sh * (sh - sqmz) / (
        fem::pow2((sh - sqmz)) + fem::pow2(gzmz));
      zz = 1.f / fem::pow2((16.f * xw * (1.f - xw))) * fem::pow2(
        sh) / (fem::pow2((sh - sqmz)) + fem::pow2(gzmz));
      if (mstp(43) == 1) {
        //C...Only gamma* production included
        gz = 0.f;
        zz = 0.f;
      }
      else if (mstp(43) == 2) {
        //C...Only Z0 production included
        gg = 0.f;
        gz = 0.f;
      }
      asym = 2.f * (ei * ai * gz * ef * af + 4.f * vi * ai * zz *
        vf * af) / (fem::pow2(ei) * gg * fem::pow2(ef) + ei * vi *
        gz * ef * vf + (fem::pow2(vi) + fem::pow2(ai)) * zz * (
        fem::pow2(vf) + fem::pow2(af)));
      wt = 1.f + asym * cthe(jt) + fem::pow2(cthe(jt));
      wtmax = 2.f + fem::abs(asym);
    }
    else {
      //C...Angular weight for gamma*/Z0 -> H+ + H-
      wt = 1.f - fem::pow2(cthe(jt));
      wtmax = 1.f;
    }
    //C
  }
  else if (isub == 2) {
    //C...Angular weight for W+/- -> 2 quarks/leptons
    wt = fem::pow2((1.f + cthe(jt)));
    wtmax = 4.f;
    //C
  }
  else if (isub == 15 || isub == 19) {
    //C...Angular weight for f + fb -> gluon/gamma + Z0 ->
    //C...-> gluon/gamma + 2 quarks/leptons
    wt = (fem::pow2((coup(1, 3) * coup(3, 3))) + fem::pow2((coup(1,
      4) * coup(3, 4)))) * (fem::pow2(pkk(1, 3)) + fem::pow2(pkk(2,
      4))) + (fem::pow2((coup(1, 3) * coup(3, 4))) + fem::pow2((coup(1,
      4) * coup(3, 3)))) * (fem::pow2(pkk(1, 4)) + fem::pow2(pkk(2,
      3)));
    wtmax = (fem::pow2(coup(1, 3)) + fem::pow2(coup(1, 4))) * (
      fem::pow2(coup(3, 3)) + fem::pow2(coup(3, 4))) * (fem::pow2((pkk(1,
      3) + pkk(1, 4))) + fem::pow2((pkk(2, 3) + pkk(2, 4))));
    //C
  }
  else if (isub == 16 || isub == 20) {
    //C...Angular weight for f + fb' -> gluon/gamma + W+/- ->
    //C...-> gluon/gamma + 2 quarks/leptons
    wt = fem::pow2(pkk(1, 3)) + fem::pow2(pkk(2, 4));
    wtmax = fem::pow2((pkk(1, 3) + pkk(1, 4))) + fem::pow2((pkk(2,
      3) + pkk(2, 4)));
    //C
  }
  else if (isub == 22) {
    //C...Angular weight for f + fb -> Z0 + Z0 -> 4 quarks/leptons
    s34 = fem::pow2(p(iref(ip, iord), 5));
    s56 = fem::pow2(p(iref(ip, 3 - iord), 5));
    ti = pkk(1, 3) + pkk(1, 4) + s34;
    ui = pkk(1, 5) + pkk(1, 6) + s56;
    wt = fem::pow4(coup(1, 3)) * (fem::pow2((coup(3, 3) * coup(5,
      3) * fem::abs(fgk(1, 2, 3, 4, 5, 6) / ti + fgk(1, 2, 5, 6, 3,
      4) / ui))) + fem::pow2((coup(3, 4) * coup(5, 3) * fem::abs(fgk(1,
      2, 4, 3, 5, 6) / ti + fgk(1, 2, 5, 6, 4, 3) / ui))) + fem::pow2((coup(3,
      3) * coup(5, 4) * fem::abs(fgk(1, 2, 3, 4, 6, 5) / ti + fgk(1,
      2, 6, 5, 3, 4) / ui))) + fem::pow2((coup(3, 4) * coup(5, 4) *
      fem::abs(fgk(1, 2, 4, 3, 6, 5) / ti + fgk(1, 2, 6, 5, 4, 3) /
      ui)))) + fem::pow4(coup(1, 4)) * (fem::pow2((coup(3, 3) * coup(5,
      3) * fem::abs(fgk(2, 1, 5, 6, 3, 4) / ti + fgk(2, 1, 3, 4, 5,
      6) / ui))) + fem::pow2((coup(3, 4) * coup(5, 3) * fem::abs(fgk(2,
      1, 6, 5, 3, 4) / ti + fgk(2, 1, 3, 4, 6, 5) / ui))) + fem::pow2((coup(3,
      3) * coup(5, 4) * fem::abs(fgk(2, 1, 5, 6, 4, 3) / ti + fgk(2,
      1, 4, 3, 5, 6) / ui))) + fem::pow2((coup(3, 4) * coup(5, 4) *
      fem::abs(fgk(2, 1, 6, 5, 4, 3) / ti + fgk(2, 1, 4, 3, 6, 5) /
      ui))));
    wtmax = 4.f * s34 * s56 * (fem::pow4(coup(1, 3)) + fem::pow4(coup(1,
      4))) * (fem::pow2(coup(3, 3)) + fem::pow2(coup(3, 4))) * (
      fem::pow2(coup(5, 3)) + fem::pow2(coup(5, 4))) * 4.f * (ti /
      ui + ui / ti + 2.f * sh * (s34 + s56) / (ti * ui) - s34 * s56 *
      (1.f / fem::pow2(ti) + 1.f / fem::pow2(ui)));
    //C
  }
  else if (isub == 23) {
    //C...Angular weight for f + fb' -> Z0 + W +/- -> 4 quarks/leptons
    d34 = fem::dble(fem::pow2(p(iref(ip, iord), 5)));
    d56 = fem::dble(fem::pow2(p(iref(ip, 3 - iord), 5)));
    dt = fem::dble(pkk(1, 3) + pkk(1, 4)) + d34;
    du = fem::dble(pkk(1, 5) + pkk(1, 6)) + d56;
    cawz = coup(2, 3) / fem::sngl(dt) - 2.f * (1.f - xw) * coup(1,
      2) / (sh - sqmw);
    cbwz = coup(1, 3) / fem::sngl(du) + 2.f * (1.f - xw) * coup(1,
      2) / (sh - sqmw);
    wt = fem::pow2(coup(5, 3)) * fem::pow2(fem::abs(cawz * fgk(1, 2, 3,
      4, 5, 6) + cbwz * fgk(1, 2, 5, 6, 3, 4))) + fem::pow2(coup(5,
      4)) * fem::pow2(fem::abs(cawz * fgk(1, 2, 3, 4, 6, 5) + cbwz * fgk(1,
      2, 6, 5, 3, 4)));
    wtmax = 4.f * fem::sngl(d34 * d56) * (fem::pow2(coup(5, 3)) +
      fem::pow2(coup(5, 4))) * (fem::pow2(cawz) * fem::sngl(digk(dt,
      du)) + fem::pow2(cbwz) * fem::sngl(digk(du, dt)) + cawz * cbwz *
      fem::sngl(djgk(dt, du)));
    //C
  }
  else if (isub == 24) {
    //C...Angular weight for f + fb -> Z0 + H0 -> 2 quarks/leptons + H0
    wt = (fem::pow2((coup(1, 3) * coup(3, 3))) + fem::pow2((coup(1,
      4) * coup(3, 4)))) * pkk(1, 3) * pkk(2, 4) + (fem::pow2((coup(1,
      3) * coup(3, 4))) + fem::pow2((coup(1, 4) * coup(3, 3)))) * pkk(1,
      4) * pkk(2, 3);
    wtmax = (fem::pow2(coup(1, 3)) + fem::pow2(coup(1, 4))) * (
      fem::pow2(coup(3, 3)) + fem::pow2(coup(3, 4))) * (pkk(1, 3) + pkk(1,
      4)) * (pkk(2, 3) + pkk(2, 4));
    //C
  }
  else if (isub == 25) {
    //C...Angular weight for f + fb -> W+ + W- -> 4 quarks/leptons
    d34 = fem::dble(fem::pow2(p(iref(ip, iord), 5)));
    d56 = fem::dble(fem::pow2(p(iref(ip, 3 - iord), 5)));
    dt = fem::dble(pkk(1, 3) + pkk(1, 4)) + d34;
    du = fem::dble(pkk(1, 5) + pkk(1, 6)) + d56;
    cdww = (coup(1, 3) * sqmz / (sh - sqmz) + coup(1, 2)) / sh;
    caww = cdww + 0.5f * (coup(1, 2) + 1.f) / fem::sngl(dt);
    cbww = cdww + 0.5f * (coup(1, 2) - 1.f) / fem::sngl(du);
    ccww = coup(1, 4) * sqmz / (sh - sqmz) / sh;
    wt = fem::pow2(fem::abs(caww * fgk(1, 2, 3, 4, 5, 6) - cbww * fgk(1,
      2, 5, 6, 3, 4))) + fem::pow2(ccww) * fem::pow2(fem::abs(fgk(2,
      1, 5, 6, 3, 4) - fgk(2, 1, 3, 4, 5, 6)));
    wtmax = 4.f * fem::sngl(d34 * d56) * (fem::pow2(caww) * fem::sngl(digk(dt,
      du)) + fem::pow2(cbww) * fem::sngl(digk(du, dt)) - caww *
      cbww * fem::sngl(djgk(dt, du)) + fem::pow2(ccww) * fem::sngl(digk(dt,
      du) + digk(du, dt) - djgk(dt, du)));
    //C
  }
  else if (isub == 26) {
    //C...Angular weight for f + fb' -> W+/- + H0 -> 2 quarks/leptons + H0
    wt = pkk(1, 3) * pkk(2, 4);
    wtmax = (pkk(1, 3) + pkk(1, 4)) * (pkk(2, 3) + pkk(2, 4));
    //C
  }
  else if (isub == 30) {
    //C...Angular weight for f + g -> f + Z0 -> f + 2 quarks/leptons
    if (k(ilin(1), 2) > 0) {
      wt = (fem::pow2((coup(1, 3) * coup(3, 3))) + fem::pow2((coup(1,
        4) * coup(3, 4)))) * (fem::pow2(pkk(1, 4)) + fem::pow2(pkk(3,
        5))) + (fem::pow2((coup(1, 3) * coup(3, 4))) + fem::pow2((coup(1,
        4) * coup(3, 3)))) * (fem::pow2(pkk(1, 3)) + fem::pow2(pkk(4,
        5)));
    }
    if (k(ilin(1), 2) < 0) {
      wt = (fem::pow2((coup(1, 3) * coup(3, 3))) + fem::pow2((coup(1,
        4) * coup(3, 4)))) * (fem::pow2(pkk(1, 3)) + fem::pow2(pkk(4,
        5))) + (fem::pow2((coup(1, 3) * coup(3, 4))) + fem::pow2((coup(1,
        4) * coup(3, 3)))) * (fem::pow2(pkk(1, 4)) + fem::pow2(pkk(3,
        5)));
    }
    wtmax = (fem::pow2(coup(1, 3)) + fem::pow2(coup(1, 4))) * (
      fem::pow2(coup(3, 3)) + fem::pow2(coup(3, 4))) * (fem::pow2((pkk(1,
      3) + pkk(1, 4))) + fem::pow2((pkk(3, 5) + pkk(4, 5))));
    //C
  }
  else if (isub == 31) {
    //C...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons
    if (k(ilin(1), 2) > 0) {
      wt = fem::pow2(pkk(1, 4)) + fem::pow2(pkk(3, 5));
    }
    if (k(ilin(1), 2) < 0) {
      wt = fem::pow2(pkk(1, 3)) + fem::pow2(pkk(4, 5));
    }
    wtmax = fem::pow2((pkk(1, 3) + pkk(1, 4))) + fem::pow2((pkk(3,
      5) + pkk(4, 5)));
    //C
  }
  else if (isub == 141) {
    //C...Angular weight for gamma*/Z0/Z'0 -> 2 quarks/leptons
    ei = kchg(fem::iabs(mint(15)), 1) / 3.f;
    ai = fem::sign(1.f, ei + 0.1f);
    vi = ai - 4.f * ei * xw;
    api = fem::sign(1.f, ei + 0.1f);
    vpi = api - 4.f * ei * xw;
    ef = kchg(kfa, 1) / 3.f;
    af = fem::sign(1.f, ef + 0.1f);
    vf = af - 4.f * ef * xw;
    apf = fem::sign(1.f, ef + 0.1f);
    vpf = apf - 4.f * ef * xw;
    gg = 1.f;
    gz = 1.f / (8.f * xw * (1.f - xw)) * sh * (sh - sqmz) / (
      fem::pow2((sh - sqmz)) + fem::pow2(gzmz));
    gzp = 1.f / (8.f * xw * (1.f - xw)) * sh * (sh - sqmzp) / (
      fem::pow2((sh - sqmzp)) + fem::pow2(gzmzp));
    zz = 1.f / fem::pow2((16.f * xw * (1.f - xw))) * fem::pow2(sh) / (
      fem::pow2((sh - sqmz)) + fem::pow2(gzmz));
    zzp = 2.f / fem::pow2((16.f * xw * (1.f - xw))) * fem::pow2(sh) *
      ((sh - sqmz) * (sh - sqmzp) + gzmz * gzmzp) / ((fem::pow2((sh -
      sqmz)) + fem::pow2(gzmz)) * (fem::pow2((sh - sqmzp)) +
      fem::pow2(gzmzp)));
    zpzp = 1.f / fem::pow2((16.f * xw * (1.f - xw))) * fem::pow2(
      sh) / (fem::pow2((sh - sqmzp)) + fem::pow2(gzmzp));
    if (mstp(44) == 1) {
      //C...Only gamma* production included
      gz = 0.f;
      gzp = 0.f;
      zz = 0.f;
      zzp = 0.f;
      zpzp = 0.f;
    }
    else if (mstp(44) == 2) {
      //C...Only Z0 production included
      gg = 0.f;
      gz = 0.f;
      gzp = 0.f;
      zzp = 0.f;
      zpzp = 0.f;
    }
    else if (mstp(44) == 3) {
      //C...Only Z'0 production included
      gg = 0.f;
      gz = 0.f;
      gzp = 0.f;
      zz = 0.f;
      zzp = 0.f;
    }
    else if (mstp(44) == 4) {
      //C...Only gamma*/Z0 production included
      gzp = 0.f;
      zzp = 0.f;
      zpzp = 0.f;
    }
    else if (mstp(44) == 5) {
      //C...Only gamma*/Z'0 production included
      gz = 0.f;
      zz = 0.f;
      zzp = 0.f;
    }
    else if (mstp(44) == 6) {
      //C...Only Z0/Z'0 production included
      gg = 0.f;
      gz = 0.f;
      gzp = 0.f;
    }
    asym = 2.f * (ei * ai * gz * ef * af + ei * api * gzp * ef *
      apf + 4.f * vi * ai * zz * vf * af + (vi * api + vpi * ai) *
      zzp * (vf * apf + vpf * af) + 4.f * vpi * api * zpzp * vpf *
      apf) / (fem::pow2(ei) * gg * fem::pow2(ef) + ei * vi * gz *
      ef * vf + ei * vpi * gzp * ef * vpf + (fem::pow2(vi) +
      fem::pow2(ai)) * zz * (fem::pow2(vf) + fem::pow2(af)) + (vi *
      vpi + ai * api) * zzp * (vf * vpf + af * apf) + (fem::pow2(
      vpi) + fem::pow2(api)) * zpzp * (fem::pow2(vpf) + fem::pow2(
      apf)));
    wt = 1.f + asym * cthe(jt) + fem::pow2(cthe(jt));
    wtmax = 2.f + fem::abs(asym);
    //C
  }
  else {
    wt = 1.f;
    wtmax = 1.f;
  }
  //C...Obtain correct angular distribution by rejection techniques.
  if (wt < rlu(cmn, 0) * wtmax) {
    goto statement_420;
  }
  //C
  //C...Construct massive four-vectors using angles chosen. Mark decayed
  //C...resonances, add documentation lines. Shower evolution.
  statement_500:
  FEM_DO_SAFE(jt, 1, jtmax) {
    if (kdcy(jt) == 0) {
      goto statement_520;
    }
    id = iref(ip, jt);
    ludbrb(nsd(jt) + 1, nsd(jt) + 2, fem::acos(cthe(jt)), phi(jt),
      fem::dble(p(id, 1) / p(id, 4)), fem::dble(p(id, 2) / p(id, 4)),
      fem::dble(p(id, 3) / p(id, 4)));
    k(id, 1) += 10;
    k(id, 4) = nsd(jt) + 1;
    k(id, 5) = nsd(jt) + 2;
    idoc = mint(83) + mint(4);
    FEM_DO_SAFE(i, nsd(jt) + 1, nsd(jt) + 2) {
      mint(4)++;
      i1 = mint(83) + mint(4);
      k(i, 3) = i1;
      k(i1, 1) = 21;
      k(i1, 2) = k(i, 2);
      k(i1, 3) = iref(ip, jt + 2);
      FEM_DO_SAFE(j, 1, 5) {
        p(i1, j) = p(i, j);
      }
    }
    if (jtmax == 1) {
      mint(7) = mint(83) + 6 + 2 * iset(isub);
      mint(8) = mint(83) + 7 + 2 * iset(isub);
    }
    //Clin-8/19/02 avoid actual argument in common blocks of LUSHOW:
    //C      IF(MSTP(71).GE.1.AND.KDCY(JT).EQ.1) CALL LUSHOW(NSD(JT)+1,
    //C     &NSD(JT)+2,P(ID,5))
    pid5 = p(id, 5);
    if (mstp(71) >= 1 && kdcy(jt) == 1) {
      lushow(cmn, nsd(jt) + 1, nsd(jt) + 2, pid5);
    }
    //C
    //C...Check if new resonances were produced, loop back if needed.
    if (kdcy(jt) != 3) {
      goto statement_520;
    }
    np++;
    iref(np, 1) = nsd(jt) + 1;
    iref(np, 2) = nsd(jt) + 2;
    iref(np, 3) = idoc + 1;
    iref(np, 4) = idoc + 2;
    iref(np, 5) = k(iref(ip, jt), 2);
    iref(np, 6) = iref(ip, jt);
    statement_520:;
  }
  statement_530:
  if (ip < np) {
    goto statement_100;
  }
  //C
}

//C
//C*********************************************************************
//C
void
pydiff(
  common& cmn)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  arr_cref<float> parj(cmn.parj, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_cref<float> vint(cmn.vint, dimension(400));
  //
  int jt = fem::int0;
  int i = fem::int0;
  int j = fem::int0;
  int isub = fem::int0;
  float sqlam = fem::float0;
  float pz = fem::float0;
  float pe = fem::float0;
  int imb = fem::int0;
  float chik = fem::float0;
  float chi = fem::float0;
  float cut = fem::float0;
  float cutr = fem::float0;
  float chir = fem::float0;
  float sqm = fem::float0;
  float pzi = fem::float0;
  float pei = fem::float0;
  float pqqp = fem::float0;
  //C
  //C...Handles diffractive and elastic scattering.
  //C
  //C...Reset K, P and V vectors. Store incoming particles.
  FEM_DO_SAFE(jt, 1, mstp(126) + 10) {
    i = mint(83) + jt;
    FEM_DO_SAFE(j, 1, 5) {
      k(i, j) = 0;
      p(i, j) = 0.f;
      v(i, j) = 0.f;
    }
  }
  n = mint(84);
  mint(3) = 0;
  mint(21) = 0;
  mint(22) = 0;
  mint(23) = 0;
  mint(24) = 0;
  mint(4) = 4;
  FEM_DO_SAFE(jt, 1, 2) {
    i = mint(83) + jt;
    k(i, 1) = 21;
    k(i, 2) = mint(10 + jt);
    p(i, 5) = vint(2 + jt);
    p(i, 3) = vint(5) * fem::pow((-1), (jt + 1));
    p(i, 4) = fem::sqrt(fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
  }
  mint(6) = 2;
  //C
  //C...Subprocess; kinematics.
  isub = mint(1);
  sqlam = fem::pow2((vint(2) - vint(63) - vint(64))) - 4.f * vint(
    63) * vint(64);
  pz = fem::sqrt(sqlam) / (2.f * vint(1));
  FEM_DO_SAFE(jt, 1, 2) {
    i = mint(83) + jt;
    pe = (vint(2) + vint(62 + jt) - vint(65 - jt)) / (2.f * vint(1));
    //C
    //C...Elastically scattered particle.
    if (mint(16 + jt) <= 0) {
      n++;
      k(n, 1) = 1;
      k(n, 2) = k(i, 2);
      k(n, 3) = i + 2;
      p(n, 3) = pz * fem::pow((-1), (jt + 1));
      p(n, 4) = pe;
      p(n, 5) = p(i, 5);
      //C
      //C...Diffracted particle: valence quark kicked out.
    }
    else if (mstp(101) == 1) {
      n += 2;
      k(n - 1, 1) = 2;
      k(n, 1) = 1;
      k(n - 1, 3) = i + 2;
      k(n, 3) = i + 2;
      pyspli(cmn, k(i, 2), 21, k(n, 2), k(n - 1, 2));
      p(n - 1, 5) = ulmass(cmn, k(n - 1, 2));
      p(n, 5) = ulmass(cmn, k(n, 2));
      sqlam = fem::pow2((vint(62 + jt) - fem::pow2(p(n - 1, 5)) -
        fem::pow2(p(n, 5)))) - 4.f * fem::pow2(p(n - 1, 5)) * fem::pow2(p(n,
        5));
      p(n - 1, 3) = (pe * fem::sqrt(sqlam) + pz * (vint(62 + jt) +
        fem::pow2(p(n - 1, 5)) - fem::pow2(p(n, 5)))) / (2.f * vint(
        62 + jt)) * fem::pow((-1), (jt + 1));
      p(n - 1, 4) = fem::sqrt(fem::pow2(p(n - 1, 3)) + fem::pow2(p(n - 1, 5)));
      p(n, 3) = pz * fem::pow((-1), (jt + 1)) - p(n - 1, 3);
      p(n, 4) = fem::sqrt(fem::pow2(p(n, 3)) + fem::pow2(p(n, 5)));
      //C
      //C...Diffracted particle: gluon kicked out.
    }
    else {
      n += 3;
      k(n - 2, 1) = 2;
      k(n - 1, 1) = 2;
      k(n, 1) = 1;
      k(n - 2, 3) = i + 2;
      k(n - 1, 3) = i + 2;
      k(n, 3) = i + 2;
      pyspli(cmn, k(i, 2), 21, k(n, 2), k(n - 2, 2));
      k(n - 1, 2) = 21;
      p(n - 2, 5) = ulmass(cmn, k(n - 2, 2));
      p(n - 1, 5) = 0.f;
      p(n, 5) = ulmass(cmn, k(n, 2));
      //C...Energy distribution for particle into two jets.
      statement_120:
      imb = 1;
      if (fem::mod(k(i, 2) / 1000, 10) != 0) {
        imb = 2;
      }
      chik = parp(92 + 2 * imb);
      if (mstp(92) <= 1) {
        if (imb == 1) {
          chi = rlu(cmn, 0);
        }
        if (imb == 2) {
          chi = 1.f - fem::sqrt(rlu(cmn, 0));
        }
      }
      else if (mstp(92) == 2) {
        chi = 1.f - fem::pow(rlu(cmn, 0), (1.f / (1.f + chik)));
      }
      else if (mstp(92) == 3) {
        cut = 2.f * 0.3f / vint(1);
        statement_130:
        chi = fem::pow2(rlu(cmn, 0));
        if (fem::pow((fem::pow2(chi) / (fem::pow2(chi) + fem::pow2(cut))),
            0.25f) * fem::pow((1.f - chi), chik) < rlu(cmn, 0)) {
          goto statement_130;
        }
      }
      else {
        cut = 2.f * 0.3f / vint(1);
        cutr = (1.f + fem::sqrt(1.f + fem::pow2(cut))) / cut;
        statement_140:
        chir = cut * fem::pow(cutr, rlu(cmn, 0));
        chi = (fem::pow2(chir) - fem::pow2(cut)) / (2.f * chir);
        if (fem::pow((1.f - chi), chik) < rlu(cmn, 0)) {
          goto statement_140;
        }
      }
      if (chi < fem::pow2(p(n, 5)) / vint(62 + jt) || chi > 1.f -
          fem::pow2(p(n - 2, 5)) / vint(62 + jt)) {
        goto statement_120;
      }
      sqm = fem::pow2(p(n - 2, 5)) / (1.f - chi) + fem::pow2(p(n, 5)) / chi;
      if (fem::pow2((fem::sqrt(sqm) + parj(32))) >= vint(62 + jt)) {
        goto statement_120;
      }
      pzi = (pe * (vint(62 + jt) - sqm) + pz * (vint(62 + jt) +
        sqm)) / (2.f * vint(62 + jt));
      pei = fem::sqrt(fem::pow2(pzi) + sqm);
      pqqp = (1.f - chi) * (pei + pzi);
      p(n - 2, 3) = 0.5f * (pqqp - fem::pow2(p(n - 2, 5)) / pqqp) *
        fem::pow((-1), (jt + 1));
      p(n - 2, 4) = fem::sqrt(fem::pow2(p(n - 2, 3)) + fem::pow2(p(n - 2, 5)));
      p(n - 1, 3) = (pz - pzi) * fem::pow((-1), (jt + 1));
      p(n - 1, 4) = fem::abs(p(n - 1, 3));
      p(n, 3) = pzi * fem::pow((-1), (jt + 1)) - p(n - 2, 3);
      p(n, 4) = fem::sqrt(fem::pow2(p(n, 3)) + fem::pow2(p(n, 5)));
    }
    //C
    //C...Documentation lines.
    k(i + 2, 1) = 21;
    if (mint(16 + jt) == 0) {
      k(i + 2, 2) = mint(10 + jt);
    }
    if (mint(16 + jt) != 0) {
      k(i + 2, 2) = 10 * (mint(10 + jt) / 10);
    }
    k(i + 2, 3) = i;
    p(i + 2, 3) = pz * fem::pow((-1), (jt + 1));
    p(i + 2, 4) = pe;
    p(i + 2, 5) = fem::sqrt(vint(62 + jt));
  }
  //C
  //C...Rotate outgoing partons/particles using cos(theta).
  ludbrb(mint(83) + 3, n, fem::acos(vint(23)), vint(24), 0e0, 0e0, 0e0);
  //C
}

//C
//C*********************************************************************
//C
void
pyfram(
  common& cmn,
  int const& iframe)
{
  common_write write(cmn);
  // COMMON ludat1
  arr_cref<int> mstu(cmn.mstu, dimension(200));
  // COMMON pypars
  arr_ref<int> msti(cmn.msti, dimension(200));
  // COMMON pyint1
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_cref<float> vint(cmn.vint, dimension(400));
  //
  //C
  //C...Performs transformations between different coordinate frames.
  //C
  if (iframe < 1 || iframe > 2) {
    write(mstu(11),
      "(1x,'Error: illegal values in subroutine PYFRAM.',1x,"
      "'No transformation performed.',/,1x,'IFRAME =',1x,i5,'; MINT(6) =',1x,"
      "i5)"),
      iframe, mint(6);
    return;
  }
  if (iframe == mint(6)) {
    return;
  }
  //C
  if (mint(6) == 1) {
    //C...Transform from fixed target or user specified frame to
    //C...CM-frame of incoming particles.
    lurobo(cmn, 0.f, 0.f, -vint(8), -vint(9), -vint(10));
    lurobo(cmn, 0.f, -vint(7), 0.f, 0.f, 0.f);
    lurobo(cmn, -vint(6), 0.f, 0.f, 0.f, 0.f);
    mint(6) = 2;
    //C
  }
  else {
    //C...Transform from particle CM-frame to fixed target or user specified
    //C...frame.
    lurobo(cmn, vint(6), vint(7), vint(8), vint(9), vint(10));
    mint(6) = 1;
  }
  msti(6) = mint(6);
  //C
}

//C
//C*********************************************************************
//C
void
pythia(
  common& cmn)
{
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_cref<float, 2> pmas(cmn.pmas, dimension(500, 4));
  arr_cref<int> msub(cmn.msub, dimension(200));
  arr_cref<int> mstp(cmn.mstp, dimension(200));
  arr_cref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> msti(cmn.msti, dimension(200));
  arr_ref<float> pari(cmn.pari, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_cref<float> vint(cmn.vint, dimension(400));
  arr_cref<int> iset(cmn.iset, dimension(200));
  arr_ref<int, 2> ngen(cmn.ngen, dim1(0, 200).dim2(3));
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  int novl = fem::int0;
  int iovl = fem::int0;
  int isub = fem::int0;
  int j = fem::int0;
  int ipu1 = fem::int0;
  int ipu2 = fem::int0;
  int nsav1 = fem::int0;
  int nsav2 = fem::int0;
  int nsav3 = fem::int0;
  int ipu3 = fem::int0;
  int ipu4 = fem::int0;
  float qmax = fem::float0;
  int i = fem::int0;
  float pt = fem::float0;
  int ngens = fem::int0;
  float xsecs = fem::float0;
  float xmaxs = fem::float0;
  float fac = fem::float0;
  int is = fem::int0;
  float pr = fem::float0;
  //C
  //C...Administers the generation of a high-pt event via calls to a number
  //C...of subroutines; also computes cross-sections.
  //C
  //C...Loop over desired number of overlayed events (normally 1).
  mint(7) = 0;
  mint(8) = 0;
  novl = 1;
  if (mstp(131) != 0) {
    pyovly(cmn, 2);
  }
  if (mstp(131) != 0) {
    novl = mint(81);
  }
  mint(83) = 0;
  mint(84) = mstp(126);
  mstu(70) = 0;
  FEM_DO_SAFE(iovl, 1, novl) {
    if (mint(84) + 100 >= mstu(4)) {
      luerrm(cmn, 11, "(PYTHIA:) no more space in LUJETS for overlayed events");
      if (mstu(21) >= 1) {
        goto statement_200;
      }
    }
    mint(82) = iovl;
    //C
    //C...Generate variables of hard scattering.
    statement_100:
    if (iovl == 1) {
      ngen(0, 2)++;
    }
    mint(31) = 0;
    mint(51) = 0;
    pyrand(cmn);
    isub = mint(1);
    if (iovl == 1) {
      ngen(isub, 2)++;
      //C
      //C...Store information on hard interaction.
      FEM_DO_SAFE(j, 1, 200) {
        msti(j) = 0;
        pari(j) = 0.f;
      }
      msti(1) = mint(1);
      msti(2) = mint(2);
      msti(11) = mint(11);
      msti(12) = mint(12);
      msti(15) = mint(15);
      msti(16) = mint(16);
      msti(17) = mint(17);
      msti(18) = mint(18);
      pari(11) = vint(1);
      pari(12) = vint(2);
      if (isub != 95) {
        FEM_DO_SAFE(j, 13, 22) {
          pari(j) = vint(30 + j);
        }
        pari(33) = vint(41);
        pari(34) = vint(42);
        pari(35) = pari(33) - pari(34);
        pari(36) = vint(21);
        pari(37) = vint(22);
        pari(38) = vint(26);
        pari(41) = vint(23);
      }
    }
    //C
    if (mstp(111) ==  - 1) {
      goto statement_160;
    }
    if (isub <= 90 || isub >= 95) {
      //C...Hard scattering (including low-pT):
      //C...reconstruct kinematics and colour flow of hard scattering.
      pyscat(cmn);
      if (mint(51) == 1) {
        goto statement_100;
      }
      //C
      //C...Showering of initial state partons (optional).
      ipu1 = mint(84) + 1;
      ipu2 = mint(84) + 2;
      if (mstp(61) >= 1 && mint(43) != 1 && isub != 95) {
        pysspa(cmn, ipu1, ipu2);
      }
      nsav1 = n;
      //C
      //C...Multiple interactions.
      if (mstp(81) >= 1 && mint(43) == 4 && isub != 95) {
        pymult(cmn, 6);
      }
      mint(1) = isub;
      nsav2 = n;
      //C
      //C...Hadron remnants and primordial kT.
      pyremn(cmn, ipu1, ipu2);
      if (mint(51) == 1) {
        goto statement_100;
      }
      nsav3 = n;
      //C
      //C...Showering of final state partons (optional).
      ipu3 = mint(84) + 3;
      ipu4 = mint(84) + 4;
      if (mstp(71) >= 1 && isub != 95 && k(ipu3, 1) > 0 && k(ipu3,
          1) <= 10 && k(ipu4, 1) > 0 && k(ipu4, 1) <= 10) {
        qmax = fem::sqrt(parp(71) * vint(52));
        if (isub == 5) {
          qmax = fem::sqrt(fem::pow2(pmas(23, 1)));
        }
        if (isub == 8) {
          qmax = fem::sqrt(fem::pow2(pmas(24, 1)));
        }
        lushow(cmn, ipu3, ipu4, qmax);
      }
      //C
      //C...Sum up transverse and longitudinal momenta.
      if (iovl == 1) {
        pari(65) = 2.f * pari(17);
        FEM_DO_SAFE(i, mstp(126) + 1, n) {
          if (k(i, 1) <= 0 || k(i, 1) > 10) {
            goto statement_130;
          }
          pt = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)));
          pari(69) += pt;
          if (i <= nsav1 || i > nsav3) {
            pari(66) += pt;
          }
          if (i > nsav1 && i <= nsav2) {
            pari(68) += pt;
          }
          statement_130:;
        }
        pari(67) = pari(68);
        pari(71) = vint(151);
        pari(72) = vint(152);
        pari(73) = vint(151);
        pari(74) = vint(152);
      }
      //C
      //C...Decay of final state resonances.
      if (mstp(41) >= 1 && isub != 95) {
        pyresd(cmn);
      }
      //C
    }
    else {
      //C...Diffractive and elastic scattering.
      pydiff(cmn);
      if (iovl == 1) {
        pari(65) = 2.f * pari(17);
        pari(66) = pari(65);
        pari(69) = pari(65);
      }
    }
    //C
    //C...Recalculate energies from momenta and masses (if desired).
    if (mstp(113) >= 1) {
      FEM_DO_SAFE(i, mint(83) + 1, n) {
        if (k(i, 1) > 0 && k(i, 1) <= 10) {
          p(i, 4) = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i,
            2)) + fem::pow2(p(i, 3)) + fem::pow2(p(i, 5)));
        }
      }
    }
    //C
    //C...Rearrange partons along strings, check invariant mass cuts.
    mstu(28) = 0;
    luprep(cmn, mint(84) + 1);
    if (mstp(112) == 1 && mstu(28) == 3) {
      goto statement_100;
    }
    if (mstp(125) == 0 || mstp(125) == 1) {
      FEM_DO_SAFE(i, mint(84) + 1, n) {
        if (k(i, 2) != 94) {
          goto statement_150;
        }
        k(i + 1, 3) = fem::mod(k(i + 1, 4) / mstu(5), mstu(5));
        k(i + 2, 3) = fem::mod(k(i + 2, 4) / mstu(5), mstu(5));
        statement_150:;
      }
      luedit(cmn, 12);
      luedit(cmn, 14);
      if (mstp(125) == 0) {
        luedit(cmn, 15);
      }
      if (mstp(125) == 0) {
        mint(4) = 0;
      }
    }
    //C
    //C...Introduce separators between sections in LULIST event listing.
    if (iovl == 1 && mstp(125) <= 0) {
      mstu(70) = 1;
      mstu(71) = n;
    }
    else if (iovl == 1) {
      mstu(70) = 3;
      mstu(71) = 2;
      mstu(72) = mint(4);
      mstu(73) = n;
    }
    //C
    //C...Perform hadronization (if desired).
    if (mstp(111) >= 1) {
      luexec(cmn);
    }
    if (mstp(125) == 0 || mstp(125) == 1) {
      luedit(cmn, 14);
    }
    //C
    //C...Calculate Monte Carlo estimates of cross-sections.
    statement_160:
    if (iovl == 1) {
      if (mstp(111) !=  - 1) {
        ngen(isub, 3)++;
      }
      ngen(0, 3)++;
      xsec(0, 3) = 0.f;
      FEM_DO_SAFE(i, 1, 200) {
        if (i == 96) {
          xsec(i, 3) = 0.f;
        }
        else if (msub(95) == 1 && (i == 11 || i == 12 || i == 13 ||
          i == 28 || i == 53 || i == 68)) {
          xsec(i, 3) = xsec(96, 2) * ngen(i, 3) / fem::max(1.f,
            fem::ffloat(ngen(96, 1)) * fem::ffloat(ngen(96, 2)));
        }
        else if (ngen(i, 1) == 0) {
          xsec(i, 3) = 0.f;
        }
        else if (ngen(i, 2) == 0) {
          xsec(i, 3) = xsec(i, 2) * ngen(0, 3) / (fem::ffloat(ngen(i,
            1)) * fem::ffloat(ngen(0, 2)));
        }
        else {
          xsec(i, 3) = xsec(i, 2) * ngen(i, 3) / (fem::ffloat(ngen(i,
            1)) * fem::ffloat(ngen(i, 2)));
        }
        xsec(0, 3) += xsec(i, 3);
      }
      if (msub(95) == 1) {
        ngens = ngen(91, 3) + ngen(92, 3) + ngen(93, 3) + ngen(94,
          3) + ngen(95, 3);
        xsecs = xsec(91, 3) + xsec(92, 3) + xsec(93, 3) + xsec(94,
          3) + xsec(95, 3);
        xmaxs = xsec(95, 1);
        if (msub(91) == 1) {
          xmaxs += xsec(91, 1);
        }
        if (msub(92) == 1) {
          xmaxs += xsec(92, 1);
        }
        if (msub(93) == 1) {
          xmaxs += xsec(93, 1);
        }
        if (msub(94) == 1) {
          xmaxs += xsec(94, 1);
        }
        fac = 1.f;
        if (ngens < ngen(0, 3)) {
          fac = (xmaxs - xsecs) / (xsec(0, 3) - xsecs);
        }
        xsec(11, 3) = fac * xsec(11, 3);
        xsec(12, 3) = fac * xsec(12, 3);
        xsec(13, 3) = fac * xsec(13, 3);
        xsec(28, 3) = fac * xsec(28, 3);
        xsec(53, 3) = fac * xsec(53, 3);
        xsec(68, 3) = fac * xsec(68, 3);
        xsec(0, 3) = xsec(91, 3) + xsec(92, 3) + xsec(93, 3) + xsec(94,
          3) + xsec(95, 1);
      }
      //C
      //C...Store final information.
      mint(5)++;
      msti(3) = mint(3);
      msti(4) = mint(4);
      msti(5) = mint(5);
      msti(6) = mint(6);
      msti(7) = mint(7);
      msti(8) = mint(8);
      msti(13) = mint(13);
      msti(14) = mint(14);
      msti(21) = mint(21);
      msti(22) = mint(22);
      msti(23) = mint(23);
      msti(24) = mint(24);
      msti(25) = mint(25);
      msti(26) = mint(26);
      msti(31) = mint(31);
      pari(1) = xsec(0, 3);
      pari(2) = xsec(0, 3) / mint(5);
      pari(31) = vint(141);
      pari(32) = vint(142);
      if (isub != 95 && mint(7) * mint(8) != 0) {
        pari(42) = 2.f * vint(47) / vint(1);
        FEM_DO_SAFE(is, 7, 8) {
          pari(36 + is) = p(mint(is), 3) / vint(1);
          pari(38 + is) = p(mint(is), 4) / vint(1);
          i = mint(is);
          pr = fem::max(1e-20f, fem::pow2(p(i, 5)) + fem::pow2(p(i,
            1)) + fem::pow2(p(i, 2)));
          pari(40 + is) = fem::sign(fem::log(fem::min((fem::sqrt(pr +
            fem::pow2(p(i, 3))) + fem::abs(p(i, 3))) / fem::sqrt(pr),
            1e20f)), p(i, 3));
          pr = fem::max(1e-20f, fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)));
          pari(42 + is) = fem::sign(fem::log(fem::min((fem::sqrt(pr +
            fem::pow2(p(i, 3))) + fem::abs(p(i, 3))) / fem::sqrt(pr),
            1e20f)), p(i, 3));
          pari(44 + is) = p(i, 3) / fem::sqrt(fem::pow2(p(i, 1)) +
            fem::pow2(p(i, 2)) + fem::pow2(p(i, 3)));
          pari(46 + is) = ulangl(cmn, p(i, 3), fem::sqrt(fem::pow2(p(i,
            1)) + fem::pow2(p(i, 2))));
          pari(48 + is) = ulangl(cmn, p(i, 1), p(i, 2));
        }
      }
      pari(61) = vint(148);
      if (iset(isub) == 1 || iset(isub) == 3) {
        mstu(161) = mint(21);
        mstu(162) = 0;
      }
      else {
        mstu(161) = mint(21);
        mstu(162) = mint(22);
      }
    }
    //C
    //C...Prepare to go to next overlayed event.
    msti(41) = iovl;
    if (iovl >= 2 && iovl <= 10) {
      msti(40 + iovl) = isub;
    }
    if (mstu(70) < 10) {
      mstu(70)++;
      mstu(70 + mstu(70)) = n;
    }
    mint(83) = n;
    mint(84) = n + mstp(126);
  }
  //C
  //C...Information on overlayed events.
  if (mstp(131) == 1 && mstp(133) >= 1) {
    pari(91) = vint(132);
    pari(92) = vint(133);
    pari(93) = vint(134);
    if (mstp(133) == 2) {
      pari(93) = pari(93) * xsec(0, 3) / vint(131);
    }
  }
  //C
  //C...Transform to the desired coordinate frame.
  statement_200:
  pyfram(cmn, mstp(124));
  //C
}

struct blockdata_pydata_save
{
};

//C
//C*********************************************************************
//C
void
blockdata_pydata(
  common& cmn)
{
  FEM_CMN_SVE(blockdata_pydata);
  // COMMON pysubs
  arr_ref<int> msub(cmn.msub, dimension(200));
  arr_ref<int, 2> kfin(cmn.kfin, dim1(2).dim2(-40, 40));
  arr_ref<float> ckin(cmn.ckin, dimension(200));
  // COMMON pypars
  arr_ref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> msti(cmn.msti, dimension(200));
  arr_ref<float> pari(cmn.pari, dimension(200));
  // COMMON pyint1
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float> vint(cmn.vint, dimension(400));
  // COMMON pyint2
  arr_ref<int> iset(cmn.iset, dimension(200));
  arr_ref<int, 2> kfpr(cmn.kfpr, dimension(200, 2));
  arr_ref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_ref<int, 3> icol(cmn.icol, dimension(40, 4, 2));
  // COMMON pyint6
  str_arr_ref<1> proc(cmn.proc, dim1(0, 200));
  //
  int i = fem::int0;
  int j = fem::int0;
  int k = fem::int0;
  if (is_called_first_time) {
    cmn.msel = 1;
    fem::data((values, 200*datum(0))), msub;
    {
      fem::data_values data((values, 40*datum(1), 0, 80*datum(1), 0,
        40*datum(1)));
      FEM_DO_SAFE(i, 1, 2) {
        FEM_DO_SAFE(j, -40, 40) {
          data, kfin(i, j);
        }
      }
    }
    {
      fem::data_values data;
      data.values, 2.0f, -1.0f, 0.0f, -1.0f, 1.0f, 1.0f, -10.f, 10.f;
      data.values, -10.f, 10.f, -10.f, 10.f, -10.f, 10.f, -10.f, 10.f;
      data.values, -1.0f, 1.0f, -1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f;
      data.values, -1.0f, 1.0f, -1.0f, 1.0f, 0.f, 0.f, 2.0f, -1.0f;
      data.values, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f;
      data.values, 160*datum(0.f);
      data, ckin;
    }
    {
      static const int values[] = {
        3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 0,
          1, 0, 3, 7, 1, 0, 0, 0, 0, 0, 1, 1, 20, 6, 0, 0, 0, 0, 0,
          0, 1, 2, 2, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 1, 1, 100, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0,
          0, 0
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 100) {
        data, mstp(i);
      }
    }
    {
      static const int values[] = {
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
          0, 1, 2, 1, 1, 20, 0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 5, 3, 1989, 11, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 101, 200) {
        data, mstp(i);
      }
    }
    {
      static const float values[] = {
        0.25f, 10.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.5f, 2.0f, 0.075f, 0.f,
          0.2f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 1.0f, 2.26f, 1.e4f, 1.e-4f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.25f, 1.0f, 0.25f, 1.0f, 2.0f,
          1.e-3f, 4.0f, 0.f, 0.f, 0.f, 4.0f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 1.6f, 1.85f, 0.5f, 0.2f, 0.33f, 0.66f,
          0.7f, 0.5f, 0.f, 0.f, 0.44f, 0.44f, 2.0f, 1.0f, 0.f, 3.0f,
          1.0f, 0.75f, 0.f, 0.f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 100) {
        data, parp(i);
      }
    }
    {
      static const float values[] = {
        -0.02f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 2.0f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.4f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.01f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
          0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 101, 200) {
        data, parp(i);
      }
    }
    fem::data((values, 200*datum(0))), msti;
    fem::data((values, 200*datum(0.f))), pari;
    fem::data((values, 400*datum(0))), mint;
    fem::data((values, 400*datum(0.f))), vint;
    {
      static const int values[] = {
        1, 1, 1, -1, 3, -1, -1, 3, -2, -2, 2, 2, 2, 2, 2, 2, -1, 2,
          2, 2, -1, 2, 2, 2, 2, 2, -1, 2, 2, 2, 2, -1, -1, -1, -1,
          -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
          -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
          -1, -1, 2, -1, -1, 4, 4, 4, -1, -1, 4, 4, -1, -1, -2, 2, 2,
          -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 0, -1, 0, 5, -2, -2,
          -2, -2
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 100) {
        data, iset(i);
      }
    }
    {
      static const int values[] = {
        -1, 1, -2, -2, -2, -2, -2, -2, -2, -2, 2, 2, 2, 2, -1, -1,
          -1, -2, -2, -2, -1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
          -2, -2, -2, -2, -2, -2, -2, -2, -2, 1, 1, 1, -2, -2, -2,
          -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 2,
          -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
          -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
          -2, -2, -2, -2, -2, -2, -2, -2, -2
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 101, 200) {
        data, iset(i);
      }
    }
    {
      static const int values[] = {
        23, 0, 24, 0, 25, 0, 24, 0, 25, 0, 24, 0, 23, 0, 25, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 21, 21, 21, 22, 21, 23, 21, 24, 21, 25,
          22, 22, 22, 23, 22, 24, 22, 25, 23, 23, 23, 24, 23, 25, 24,
          24, 24, 25, 25, 25, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0,
          21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 21, 0, 22, 0, 23, 0, 24,
          0, 25, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 21, 0, 22, 0,
          23
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 50) {
        FEM_DO_SAFE(j, 1, 2) {
          data, kfpr(i, j);
        }
      }
    }
    {
      static const int values[] = {
        0, 24, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 21, 24, 24,
          22, 24, 23, 23, 24, 24, 23, 24, 23, 25, 22, 22, 23, 23, 24,
          24, 24, 25, 25, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 51, 100) {
        FEM_DO_SAFE(j, 1, 2) {
          data, kfpr(i, j);
        }
      }
    }
    {
      static const int values[] = {
        23, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          21, 25, 0, 25, 21, 25, 22, 22, 22, 23, 23, 23, 24, 24, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 32, 0, 37, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 101, 150) {
        FEM_DO_SAFE(j, 1, 2) {
          data, kfpr(i, j);
        }
      }
    }
    {
      static const int values[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 151, 200) {
        FEM_DO_SAFE(j, 1, 2) {
          data, kfpr(i, j);
        }
      }
    }
    fem::data((values, 4000*datum(0.f))), coef;
    {
      static const int values[] = {
        4, 0, 3, 0, 2, 0, 1, 0, 3, 0, 4, 0, 1, 0, 2, 0, 2, 0, 0, 1,
          4, 0, 0, 3, 3, 0, 0, 4, 1, 0, 0, 2, 3, 0, 0, 4, 1, 4, 3, 2,
          4, 0, 0, 3, 4, 2, 1, 3, 2, 0, 4, 1, 4, 0, 2, 3, 4, 0, 3, 4,
          2, 0, 1, 2, 3, 2, 1, 0, 1, 4, 3, 0, 4, 3, 3, 0, 2, 1, 1, 0,
          3, 2, 1, 4, 1, 0, 0, 2, 2, 4, 3, 1, 2, 0, 0, 1, 3, 2, 1, 4,
          1, 4, 3, 2, 4, 2, 1, 3, 4, 2, 1, 3, 3, 4, 4, 3, 1, 2, 2, 1,
          2, 0, 3, 1, 2, 0, 0, 0, 4, 2, 1, 0, 0, 0, 1, 0, 3, 0, 0, 3,
          1, 2, 0, 0, 4, 0, 0, 4, 0, 0, 1, 2, 2, 0, 0, 1, 4, 4, 3, 3,
          2, 2, 1, 1, 4, 4, 3, 3, 3, 3, 4, 4, 1, 1, 2, 2, 3, 2, 1, 3,
          1, 2, 0, 0, 4, 2, 1, 4, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 40) {
        FEM_DO_SAFE(j, 1, 4) {
          FEM_DO_SAFE(k, 1, 2) {
            data, icol(i, j, k);
          }
        }
      }
    }
    proc(0) = "All included subprocesses   ";
    {
      static const char* values[] = {
        "f + fb -> gamma*/Z0         ",
          "f + fb' -> W+/-             ",
          "f + fb -> H0                ",
          "gamma + W+/- -> W+/-        ",
          "Z0 + Z0 -> H0               ",
          "Z0 + W+/- -> W+/-           ",
          "                            ",
          "W+ + W- -> H0               ",
          "                            ",
          "                            ",
          "f + f' -> f + f'            ",
          "f + fb -> f' + fb'          ",
          "f + fb -> g + g             ",
          "f + fb -> g + gamma         ",
          "f + fb -> g + Z0            ",
          "f + fb' -> g + W+/-         ",
          "f + fb -> g + H0            ",
          "f + fb -> gamma + gamma     ",
          "f + fb -> gamma + Z0        ",
          "f + fb' -> gamma + W+/-     "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 1, 20) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "f + fb -> gamma + H0        ",
          "f + fb -> Z0 + Z0           ",
          "f + fb' -> Z0 + W+/-        ",
          "f + fb -> Z0 + H0           ",
          "f + fb -> W+ + W-           ",
          "f + fb' -> W+/- + H0        ",
          "f + fb -> H0 + H0           ",
          "f + g -> f + g              ",
          "f + g -> f + gamma          ",
          "f + g -> f + Z0             ",
          "f + g -> f' + W+/-          ",
          "f + g -> f + H0             ",
          "f + gamma -> f + g          ",
          "f + gamma -> f + gamma      ",
          "f + gamma -> f + Z0         ",
          "f + gamma -> f' + W+/-      ",
          "f + gamma -> f + H0         ",
          "f + Z0 -> f + g             ",
          "f + Z0 -> f + gamma         ",
          "f + Z0 -> f + Z0            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 21, 40) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "f + Z0 -> f' + W+/-         ",
          "f + Z0 -> f + H0            ",
          "f + W+/- -> f' + g          ",
          "f + W+/- -> f' + gamma      ",
          "f + W+/- -> f' + Z0         ",
          "f + W+/- -> f' + W+/-       ",
          "f + W+/- -> f' + H0         ",
          "f + H0 -> f + g             ",
          "f + H0 -> f + gamma         ",
          "f + H0 -> f + Z0            ",
          "f + H0 -> f' + W+/-         ",
          "f + H0 -> f + H0            ",
          "g + g -> f + fb             ",
          "g + gamma -> f + fb         ",
          "g + Z0 -> f + fb            ",
          "g + W+/- -> f + fb'         ",
          "g + H0 -> f + fb            ",
          "gamma + gamma -> f + fb     ",
          "gamma + Z0 -> f + fb        ",
          "gamma + W+/- -> f + fb'     "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 41, 60) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "gamma + H0 -> f + fb        ",
          "Z0 + Z0 -> f + fb           ",
          "Z0 + W+/- -> f + fb'        ",
          "Z0 + H0 -> f + fb           ",
          "W+ + W- -> f + fb           ",
          "W+/- + H0 -> f + fb'        ",
          "H0 + H0 -> f + fb           ",
          "g + g -> g + g              ",
          "gamma + gamma -> W+ + W-    ",
          "gamma + W+/- -> gamma + W+/-",
          "Z0 + Z0 -> Z0 + Z0          ",
          "Z0 + Z0 -> W+ + W-          ",
          "Z0 + W+/- -> Z0 + W+/-      ",
          "Z0 + Z0 -> Z0 + H0          ",
          "W+ + W- -> gamma + gamma    ",
          "W+ + W- -> Z0 + Z0          ",
          "W+/- + W+/- -> W+/- + W+/-  ",
          "W+/- + H0 -> W+/- + H0      ",
          "H0 + H0 -> H0 + H0          ",
          "                            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 61, 80) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "q + qb -> Q + QB, massive   ",
          "g + g -> Q + QB, massive    ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "Elastic scattering          ",
          "Single diffractive          ",
          "Double diffractive          ",
          "Central diffractive         ",
          "Low-pT scattering           ",
          "Semihard QCD 2 -> 2         ",
          "                            ",
          "                            ",
          "                            ",
          "                            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 81, 100) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "g + g -> gamma*/Z0          ",
          "g + g -> H0                 ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "f + fb -> g + H0            ",
          "q + g -> q + H0             ",
          "g + g -> g + H0             ",
          "g + g -> gamma + gamma      ",
          "g + g -> gamma + Z0         ",
          "g + g -> Z0 + Z0            ",
          "g + g -> W+ + W-            ",
          "                            ",
          "                            ",
          "                            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 101, 120) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "g + g -> f + fb + H0        ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 121, 140) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "f + fb -> gamma*/Z0/Z'0     ",
          "f + fb' -> H+/-             ",
          "f + fb -> R                 ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 141, 160) {
        data, proc(i);
      }
    }
    {
      static const char* values[] = {
        "f + g -> f' + H+/-          ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            ",
          "                            "
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(i, 161, 180) {
        data, proc(i);
      }
    }
    {
      fem::data_values data((values, 20*datum("                            ")));
      FEM_DO_SAFE(i, 181, 200) {
        data, proc(i);
      }
    }
  }
  //C
  //C...Give sensible default values to all status codes and parameters.
  //C
  //C...Default values for allowed processes and kinematics constraints.
  //C
  //C...Default values for main switches and parameters. Reset information.
  //C
  //C...Constants for the generation of the various processes.
  //C
  //C...Character constants: name of processes.
  //C
}

} // namespace AMPT
