#include <fem.hpp> // Fortran EMulation library of fable module

namespace AMPT {

using namespace fem::major_types;

void
dcosh(...)
{
  throw std::runtime_error(
    "Missing function implementation: dcosh");
}

void
embedhighpt(...)
{
  throw std::runtime_error(
    "Missing function implementation: embedhighpt");
}

void
hirobo(...)
{
  throw std::runtime_error(
    "Missing function implementation: hirobo");
}

void
hjana1(...)
{
  throw std::runtime_error(
    "Missing function implementation: hjana1");
}

void
hjana2(...)
{
  throw std::runtime_error(
    "Missing function implementation: hjana2");
}

void
htop(...)
{
  throw std::runtime_error(
    "Missing function implementation: htop");
}

void
luedit(...)
{
  throw std::runtime_error(
    "Missing function implementation: luedit");
}

void
luexec(...)
{
  throw std::runtime_error(
    "Missing function implementation: luexec");
}

void
lugive(...)
{
  throw std::runtime_error(
    "Missing function implementation: lugive");
}

void
minijet_out(...)
{
  throw std::runtime_error(
    "Missing function implementation: minijet_out");
}

void
ptoh(...)
{
  throw std::runtime_error(
    "Missing function implementation: ptoh");
}

void
pyinit(...)
{
  throw std::runtime_error(
    "Missing function implementation: pyinit");
}

void
pythia(...)
{
  throw std::runtime_error(
    "Missing function implementation: pythia");
}

void
ranart(...)
{
  throw std::runtime_error(
    "Missing function implementation: ranart");
}

void
ulangl(...)
{
  throw std::runtime_error(
    "Missing function implementation: ulangl");
}

void
ulmass(...)
{
  throw std::runtime_error(
    "Missing function implementation: ulmass");
}

void
zpcmn(...)
{
  throw std::runtime_error(
    "Missing function implementation: zpcmn");
}

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

struct common_hjjet1
{
  arr<int> npj;
  arr<int, 2> kfpj;
  arr<float, 2> pjpx;
  arr<float, 2> pjpy;
  arr<float, 2> pjpz;
  arr<float, 2> pjpe;
  arr<float, 2> pjpm;
  arr<int> ntj;
  arr<int, 2> kftj;
  arr<float, 2> pjtx;
  arr<float, 2> pjty;
  arr<float, 2> pjtz;
  arr<float, 2> pjte;
  arr<float, 2> pjtm;

  common_hjjet1() :
    npj(dimension(300), fem::fill0),
    kfpj(dimension(300, 500), fem::fill0),
    pjpx(dimension(300, 500), fem::fill0),
    pjpy(dimension(300, 500), fem::fill0),
    pjpz(dimension(300, 500), fem::fill0),
    pjpe(dimension(300, 500), fem::fill0),
    pjpm(dimension(300, 500), fem::fill0),
    ntj(dimension(300), fem::fill0),
    kftj(dimension(300, 500), fem::fill0),
    pjtx(dimension(300, 500), fem::fill0),
    pjty(dimension(300, 500), fem::fill0),
    pjtz(dimension(300, 500), fem::fill0),
    pjte(dimension(300, 500), fem::fill0),
    pjtm(dimension(300, 500), fem::fill0)
  {}
};

struct common_hjjet2
{
  static const int maxstr = 150001;

  int nsg;
  arr<int> njsg;
  arr<int, 2> iasg;
  arr<int, 2> k1sg;
  arr<int, 2> k2sg;
  arr<float, 2> pxsg;
  arr<float, 2> pysg;
  arr<float, 2> pzsg;
  arr<float, 2> pesg;
  arr<float, 2> pmsg;

  common_hjjet2() :
    nsg(fem::int0),
    njsg(dimension(maxstr), fem::fill0),
    iasg(dimension(maxstr, 3), fem::fill0),
    k1sg(dimension(maxstr, 100), fem::fill0),
    k2sg(dimension(maxstr, 100), fem::fill0),
    pxsg(dimension(maxstr, 100), fem::fill0),
    pysg(dimension(maxstr, 100), fem::fill0),
    pzsg(dimension(maxstr, 100), fem::fill0),
    pesg(dimension(maxstr, 100), fem::fill0),
    pmsg(dimension(maxstr, 100), fem::fill0)
  {}
};

const int common_hjjet2::maxstr;

struct common_hstrng
{
  arr<int, 2> nfp;
  arr<float, 2> pp;
  arr<int, 2> nft;
  arr<float, 2> pt;

  common_hstrng() :
    nfp(dimension(300, 15), fem::fill0),
    pp(dimension(300, 15), fem::fill0),
    nft(dimension(300, 15), fem::fill0),
    pt(dimension(300, 15), fem::fill0)
  {}
};

struct common_rndf77
{
  int nseed;

  common_rndf77() :
    nseed(fem::int0)
  {}
};

struct common_hijdat
{
  arr<float, 2> hidat0;
  arr<float> hidat;

  common_hijdat() :
    hidat0(dimension(10, 10), fem::fill0),
    hidat(dimension(10), fem::fill0)
  {}
};

struct common_anim
{
  int nevent;
  int isoft;
  int isflag;
  int izpc;

  common_anim() :
    nevent(fem::int0),
    isoft(fem::int0),
    isflag(fem::int0),
    izpc(fem::int0)
  {}
};

struct common_hjjet4
{
  static const int maxstr = 150001;

  int ndr;
  arr<int, 2> iadr;
  arr<int> kfdr;
  arr<float, 2> pdr;

  common_hjjet4() :
    ndr(fem::int0),
    iadr(dimension(maxstr, 2), fem::fill0),
    kfdr(dimension(maxstr), fem::fill0),
    pdr(dimension(maxstr, 5), fem::fill0)
  {}
};

const int common_hjjet4::maxstr;

struct common_xydr
{
  static const int maxstr = 150001;

  arr<float, 2> rtdr;

  common_xydr() :
    rtdr(dimension(maxstr, 2), fem::fill0)
  {}
};

const int common_xydr::maxstr;

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

struct common_pyint5
{
  arr<int, 2> ngen;
  arr<float, 2> xsec;

  common_pyint5() :
    ngen(dim1(0, 200).dim2(3), fem::fill0),
    xsec(dim1(0, 200).dim2(3), fem::fill0)
  {}
};

struct common_hpint
{
  int mint4;
  int mint5;
  arr<float, 2> atco;
  arr<float> atxs;

  common_hpint() :
    mint4(fem::int0),
    mint5(fem::int0),
    atco(dimension(200, 20), fem::fill0),
    atxs(dim1(0, 200), fem::fill0)
  {}
};

struct common_phidcy
{
  int iphidcy;
  float pttrig;
  int ntrig;
  int maxmiss;
  int ipi0dcy;

  common_phidcy() :
    iphidcy(fem::int0),
    pttrig(fem::float0),
    ntrig(fem::int0),
    maxmiss(fem::int0),
    ipi0dcy(fem::int0)
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

struct common_besel
{
  float x4;

  common_besel() :
    x4(fem::float0)
  {}
};

struct common_hijhb
{
  arr<float, 2> rr;
  arr<float, 2> xx;

  common_hijhb() :
    rr(dimension(10, 201), fem::fill0),
    xx(dimension(10, 201), fem::fill0)
  {}
};

struct common_dpmcm1
{
  int jjp;
  int jjt;
  float amp;
  float amt;
  float apx0;
  float atx0;
  float ampn;
  float amtn;
  float amp0;
  float amt0;
  int nfdp;
  int nfdt;
  float wp;
  float wm;
  float sw;
  float xremp;
  float xremt;
  float dpkc1;
  float dpkc2;
  float pp11;
  float pp12;
  float pt11;
  float pt12;
  float ptp2;
  float ptt2;

  common_dpmcm1() :
    jjp(fem::int0),
    jjt(fem::int0),
    amp(fem::float0),
    amt(fem::float0),
    apx0(fem::float0),
    atx0(fem::float0),
    ampn(fem::float0),
    amtn(fem::float0),
    amp0(fem::float0),
    amt0(fem::float0),
    nfdp(fem::int0),
    nfdt(fem::int0),
    wp(fem::float0),
    wm(fem::float0),
    sw(fem::float0),
    xremp(fem::float0),
    xremt(fem::float0),
    dpkc1(fem::float0),
    dpkc2(fem::float0),
    pp11(fem::float0),
    pp12(fem::float0),
    pt11(fem::float0),
    pt12(fem::float0),
    ptp2(fem::float0),
    ptt2(fem::float0)
  {}
};

struct common_dpmcm2
{
  int ndpm;
  arr<int, 2> kdpm;
  arr<float, 2> pdpm1;
  arr<float, 2> pdpm2;

  common_dpmcm2() :
    ndpm(fem::int0),
    kdpm(dimension(20, 2), fem::fill0),
    pdpm1(dimension(20, 5), fem::fill0),
    pdpm2(dimension(20, 5), fem::fill0)
  {}
};

struct common_hjglbr
{
  int nelt;
  int ninthj;
  int nelp;
  int ninp;

  common_hjglbr() :
    nelt(fem::int0),
    ninthj(fem::int0),
    nelp(fem::int0),
    ninp(fem::int0)
  {}
};

struct common_hmain1
{
  float eatt;
  int jatt;
  int natt;
  int nt;
  int np;
  int n0;
  int n01;
  int n10;
  int n11;

  common_hmain1() :
    eatt(fem::float0),
    jatt(fem::int0),
    natt(fem::int0),
    nt(fem::int0),
    np(fem::int0),
    n0(fem::int0),
    n01(fem::int0),
    n10(fem::int0),
    n11(fem::int0)
  {}
};

struct common_hmain2
{
  static const int maxstr = 150001;

  arr<int, 2> katt;
  arr<float, 2> patt;

  common_hmain2() :
    katt(dimension(maxstr, 4), fem::fill0),
    patt(dimension(maxstr, 4), fem::fill0)
  {}
};

const int common_hmain2::maxstr;

struct common_arprc
{
  static const int maxstr = 150001;

  arr<int> itypar;
  arr<float> gxar;
  arr<float> gyar;
  arr<float> gzar;
  arr<float> ftar;
  arr<float> pxar;
  arr<float> pyar;
  arr<float> pzar;
  arr<float> pear;
  arr<float> xmar;

  common_arprc() :
    itypar(dimension(maxstr), fem::fill0),
    gxar(dimension(maxstr), fem::fill0),
    gyar(dimension(maxstr), fem::fill0),
    gzar(dimension(maxstr), fem::fill0),
    ftar(dimension(maxstr), fem::fill0),
    pxar(dimension(maxstr), fem::fill0),
    pyar(dimension(maxstr), fem::fill0),
    pzar(dimension(maxstr), fem::fill0),
    pear(dimension(maxstr), fem::fill0),
    xmar(dimension(maxstr), fem::fill0)
  {}
};

const int common_arprc::maxstr;

struct common_para1
{
  int mul;

  common_para1() :
    mul(fem::int0)
  {}
};

struct common_prec1
{
  static const int maxptn = 400001;

  arr<double> gx0;
  arr<double> gy0;
  arr<double> gz0;
  arr<double> ft0;
  arr<double> px0;
  arr<double> py0;
  arr<double> pz0;
  arr<double> e0;
  arr<double> xmass0;
  arr<int> ityp0;

  common_prec1() :
    gx0(dimension(maxptn), fem::fill0),
    gy0(dimension(maxptn), fem::fill0),
    gz0(dimension(maxptn), fem::fill0),
    ft0(dimension(maxptn), fem::fill0),
    px0(dimension(maxptn), fem::fill0),
    py0(dimension(maxptn), fem::fill0),
    pz0(dimension(maxptn), fem::fill0),
    e0(dimension(maxptn), fem::fill0),
    xmass0(dimension(maxptn), fem::fill0),
    ityp0(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec1::maxptn;

struct common_prec2
{
  static const int maxptn = 400001;

  arr<double> gx5;
  arr<double> gy5;
  arr<double> gz5;
  arr<double> ft5;
  arr<double> px5;
  arr<double> py5;
  arr<double> pz5;
  arr<double> e5;
  arr<double> xmass5;
  arr<int> ityp5;

  common_prec2() :
    gx5(dimension(maxptn), fem::fill0),
    gy5(dimension(maxptn), fem::fill0),
    gz5(dimension(maxptn), fem::fill0),
    ft5(dimension(maxptn), fem::fill0),
    px5(dimension(maxptn), fem::fill0),
    py5(dimension(maxptn), fem::fill0),
    pz5(dimension(maxptn), fem::fill0),
    e5(dimension(maxptn), fem::fill0),
    xmass5(dimension(maxptn), fem::fill0),
    ityp5(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec2::maxptn;

struct common_ilist7
{
  static const int maxptn = 400001;

  arr<int> lstrg0;
  arr<int> lpart0;

  common_ilist7() :
    lstrg0(dimension(maxptn), fem::fill0),
    lpart0(dimension(maxptn), fem::fill0)
  {}
};

const int common_ilist7::maxptn;

struct common_ilist8
{
  static const int maxptn = 400001;

  arr<int> lstrg1;
  arr<int> lpart1;

  common_ilist8() :
    lstrg1(dimension(maxptn), fem::fill0),
    lpart1(dimension(maxptn), fem::fill0)
  {}
};

const int common_ilist8::maxptn;

struct common_srec1
{
  int nsp;
  int nst;
  int nsi;

  common_srec1() :
    nsp(fem::int0),
    nst(fem::int0),
    nsi(fem::int0)
  {}
};

struct common_srec2
{
  static const int maxstr = 150001;

  arr<double> ataui;
  arr<double> zt1;
  arr<double> zt2;
  arr<double> zt3;

  common_srec2() :
    ataui(dimension(maxstr), fem::fill0),
    zt1(dimension(maxstr), fem::fill0),
    zt2(dimension(maxstr), fem::fill0),
    zt3(dimension(maxstr), fem::fill0)
  {}
};

const int common_srec2::maxstr;

struct common_frzout
{
  arr<double> xnprod;
  arr<double> etprod;
  arr<double> xnfrz;
  arr<double> etfrz;
  arr<double> dnprod;
  arr<double> detpro;
  arr<double> dnfrz;
  arr<double> detfrz;

  common_frzout() :
    xnprod(dimension(30), fem::fill0),
    etprod(dimension(30), fem::fill0),
    xnfrz(dimension(30), fem::fill0),
    etfrz(dimension(30), fem::fill0),
    dnprod(dimension(30), fem::fill0),
    detpro(dimension(30), fem::fill0),
    dnfrz(dimension(30), fem::fill0),
    detfrz(dimension(30), fem::fill0)
  {}
};

struct common_soft
{
  static const int maxstr = 150001;

  arr<double, 2> pxsgs;
  arr<double, 2> pysgs;
  arr<double, 2> pzsgs;
  arr<double, 2> pesgs;
  arr<double, 2> pmsgs;
  arr<double, 2> gxsgs;
  arr<double, 2> gysgs;
  arr<double, 2> gzsgs;
  arr<double, 2> ftsgs;
  arr<int, 2> k1sgs;
  arr<int, 2> k2sgs;
  arr<int> njsgs;

  common_soft() :
    pxsgs(dimension(maxstr, 3), fem::fill0),
    pysgs(dimension(maxstr, 3), fem::fill0),
    pzsgs(dimension(maxstr, 3), fem::fill0),
    pesgs(dimension(maxstr, 3), fem::fill0),
    pmsgs(dimension(maxstr, 3), fem::fill0),
    gxsgs(dimension(maxstr, 3), fem::fill0),
    gysgs(dimension(maxstr, 3), fem::fill0),
    gzsgs(dimension(maxstr, 3), fem::fill0),
    ftsgs(dimension(maxstr, 3), fem::fill0),
    k1sgs(dimension(maxstr, 3), fem::fill0),
    k2sgs(dimension(maxstr, 3), fem::fill0),
    njsgs(dimension(maxstr), fem::fill0)
  {}
};

const int common_soft::maxstr;

struct common_noprec
{
  static const int maxidl = 4001;

  int nnozpc;
  arr<int> itypn;
  arr<float> gxn;
  arr<float> gyn;
  arr<float> gzn;
  arr<float> ftn;
  arr<float> pxn;
  arr<float> pyn;
  arr<float> pzn;
  arr<float> een;
  arr<float> xmn;

  common_noprec() :
    nnozpc(fem::int0),
    itypn(dimension(maxidl), fem::fill0),
    gxn(dimension(maxidl), fem::fill0),
    gyn(dimension(maxidl), fem::fill0),
    gzn(dimension(maxidl), fem::fill0),
    ftn(dimension(maxidl), fem::fill0),
    pxn(dimension(maxidl), fem::fill0),
    pyn(dimension(maxidl), fem::fill0),
    pzn(dimension(maxidl), fem::fill0),
    een(dimension(maxidl), fem::fill0),
    xmn(dimension(maxidl), fem::fill0)
  {}
};

const int common_noprec::maxidl;

struct common_lastt
{
  int itimeh;
  float bimp;

  common_lastt() :
    itimeh(fem::int0),
    bimp(fem::float0)
  {}
};

struct common_arevt
{
  int iaevt;
  int iarun;
  int miss;

  common_arevt() :
    iaevt(fem::int0),
    iarun(fem::int0),
    miss(fem::int0)
  {}
};

struct common_para7
{
  int ioscar;
  int nsmbbbar;
  int nsmmeson;

  common_para7() :
    ioscar(fem::int0),
    nsmbbbar(fem::int0),
    nsmmeson(fem::int0)
  {}
};

struct common_phihj
{
  int iphirp;
  float phirp;

  common_phihj() :
    iphirp(fem::int0),
    phirp(fem::float0)
  {}
};

struct common_precpa
{
  static const int maxptn = 400001;

  arr<double> vxp0;
  arr<double> vyp0;
  arr<double> vzp0;
  arr<double> xstrg0;
  arr<double> ystrg0;
  arr<double> xstrg;
  arr<double> ystrg;
  arr<int> istrg0;
  arr<int> istrg;

  common_precpa() :
    vxp0(dimension(maxptn), fem::fill0),
    vyp0(dimension(maxptn), fem::fill0),
    vzp0(dimension(maxptn), fem::fill0),
    xstrg0(dimension(maxptn), fem::fill0),
    ystrg0(dimension(maxptn), fem::fill0),
    xstrg(dimension(maxptn), fem::fill0),
    ystrg(dimension(maxptn), fem::fill0),
    istrg0(dimension(maxptn), fem::fill0),
    istrg(dimension(maxptn), fem::fill0)
  {}
};

const int common_precpa::maxptn;

struct common_wood
{
  float r;
  float d;
  float fnorm;
  float w;

  common_wood() :
    r(fem::float0),
    d(fem::float0),
    fnorm(fem::float0),
    w(fem::float0)
  {}
};

struct common_njet
{
  int n;
  int ipcrs;

  common_njet() :
    n(fem::int0),
    ipcrs(fem::int0)
  {}
};

struct common_cmsflag
{
  double dshadow;
  int ishadow;

  common_cmsflag() :
    dshadow(fem::double0),
    ishadow(fem::int0)
  {}
};

struct common_sedvax
{
  int num1;

  common_sedvax() :
    num1(fem::int0)
  {}
};

struct common_bveg1
{
  arr<double> xl;
  arr<double> xu;
  double acc;
  int ndim;
  int ncall;
  int itmx;
  int nprn;

  common_bveg1() :
    xl(dimension(10), fem::fill0),
    xu(dimension(10), fem::fill0),
    acc(fem::double0),
    ndim(fem::int0),
    ncall(fem::int0),
    itmx(fem::int0),
    nprn(fem::int0)
  {}
};

struct common_bveg2
{
  arr<double, 2> xi;
  double si;
  double si2;
  double swgt;
  double schi;
  int ndo;
  int it;

  common_bveg2() :
    xi(dimension(50, 10), fem::fill0),
    si(fem::double0),
    si2(fem::double0),
    swgt(fem::double0),
    schi(fem::double0),
    ndo(fem::int0),
    it(fem::int0)
  {}
};

struct common_bveg3
{
  double f;
  double ti;
  double tsi;

  common_bveg3() :
    f(fem::double0),
    ti(fem::double0),
    tsi(fem::double0)
  {}
};

struct common :
  fem::common,
  common_lujets,
  common_ludat1,
  common_hparnt,
  common_hjcrdn,
  common_hjjet1,
  common_hjjet2,
  common_hstrng,
  common_rndf77,
  common_hijdat,
  common_anim,
  common_hjjet4,
  common_xydr,
  common_pysubs,
  common_pypars,
  common_pyint1,
  common_pyint2,
  common_pyint5,
  common_hpint,
  common_phidcy,
  common_ludat3,
  common_besel,
  common_hijhb,
  common_dpmcm1,
  common_dpmcm2,
  common_hjglbr,
  common_hmain1,
  common_hmain2,
  common_arprc,
  common_para1,
  common_prec1,
  common_prec2,
  common_ilist7,
  common_ilist8,
  common_srec1,
  common_srec2,
  common_frzout,
  common_soft,
  common_noprec,
  common_lastt,
  common_arevt,
  common_para7,
  common_phihj,
  common_precpa,
  common_wood,
  common_njet,
  common_cmsflag,
  common_sedvax,
  common_bveg1,
  common_bveg2,
  common_bveg3
{
  fem::cmn_sve hboost_sve;
  fem::cmn_sve quench_sve;
  fem::cmn_sve ar3jet_sve;
  fem::cmn_sve arorie_sve;
  fem::cmn_sve atrobo_sve;
  fem::cmn_sve attrad_sve;
  fem::cmn_sve hijfrg_sve;
  fem::cmn_sve hijhrd_sve;
  fem::cmn_sve jetini_sve;
  fem::cmn_sve attflv_sve;
  fem::cmn_sve hijini_sve;
  fem::cmn_sve hijels_sve;
  fem::cmn_sve gauss2_sve;
  fem::cmn_sve romg_sve;
  fem::cmn_sve hijcsc_sve;
  fem::cmn_sve hirnd2_sve;
  fem::cmn_sve hijsft_sve;
  fem::cmn_sve hirnd_sve;
  fem::cmn_sve hijing_sve;
  fem::cmn_sve ftot_sve;
  fem::cmn_sve fhin_sve;
  fem::cmn_sve ftotjt_sve;
  fem::cmn_sve ftotrg_sve;
  fem::cmn_sve sgmin_sve;
  fem::cmn_sve fnjet_sve;
  fem::cmn_sve gauss1_sve;
  fem::cmn_sve hifun_sve;
  fem::cmn_sve hijwds_sve;
  fem::cmn_sve gmre_sve;
  fem::cmn_sve parton_sve;
  fem::cmn_sve g_sve;
  fem::cmn_sve fjet_sve;
  fem::cmn_sve ghvq_sve;
  fem::cmn_sve gphotn_sve;
  fem::cmn_sve fjetrg_sve;
  fem::cmn_sve aran9_sve;
  fem::cmn_sve vegas_sve;
  fem::cmn_sve crsjet_sve;
  fem::cmn_sve hijcrs_sve;
  fem::cmn_sve hijset_sve;
  fem::cmn_sve blockdata_hidata_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct hboost_save
{
  float amt;
  double db;
  double dbeta;
  double dga;
  double dp3;
  double dp4;
  int i;
  float y;

  hboost_save() :
    amt(fem::float0),
    db(fem::double0),
    dbeta(fem::double0),
    dga(fem::double0),
    dp3(fem::double0),
    dp4(fem::double0),
    i(fem::int0),
    y(fem::float0)
  {}
};

void
hboost(
  common& cmn)
{
  FEM_CMN_SVE(hboost);
  common_write write(cmn);
  int& n = static_cast<common_lujets&>(cmn).n;
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  float& amt = sve.amt;
  double& db = sve.db;
  double& dbeta = sve.dbeta;
  double& dga = sve.dga;
  double& dp3 = sve.dp3;
  double& dp4 = sve.dp4;
  int& i = sve.i;
  float& y = sve.y;
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /LUDAT1/
  //Cc      SAVE /HPARNT/
  FEM_DO_SAFE(i, 1, n) {
    dbeta = fem::dble(p(i, 3) / p(i, 4));
    if (fem::abs(dbeta) >= 1.e0) {
      db = fem::dble(hint1(2));
      if (db > 0.99999999e0) {
        //C                ********Rescale boost vector if too close to unity.
        write(6, star), "(HIBOOT:) boost vector too large";
        db = 0.99999999e0;
      }
      dga = 1e0 / fem::sqrt(1e0 - fem::pow2(db));
      dp3 = fem::dble(p(i, 3));
      dp4 = fem::dble(p(i, 4));
      p(i, 3) = fem::sngl((dp3 + db * dp4) * dga);
      p(i, 4) = fem::sngl((dp4 + db * dp3) * dga);
      goto statement_100;
    }
    y = 0.5f * fem::sngl(fem::dlog((1.e0 + dbeta) / (1.e0 - dbeta)));
    amt = fem::sqrt(fem::pow2(p(i, 1)) + fem::pow2(p(i, 2)) + fem::pow2(p(i,
      5)));
    p(i, 3) = amt * fem::sinh(y + hint1(3));
    p(i, 4) = amt * fem::cosh(y + hint1(3));
    statement_100:;
  }
}

struct quench_save
{
  float amshu;
  float bb;
  float bbx;
  float bby;
  float de;
  float dp;
  float dp1;
  float dp2;
  float dp3;
  float dphi;
  float drr;
  float dx;
  float dy;
  float ershu;
  int i;
  int i2;
  int isg;
  int j2;
  int jp;
  int jt;
  int kp;
  int kt;
  int lq;
  arr<int> lqp;
  arr<int> lqt;
  int mp;
  int mt;
  int nq;
  float phi;
  float phip;
  float phiq;
  float phit;
  float prshu;
  float ptjet0;
  float ptot;
  float r0;
  float rd;
  float rd0;
  arr<float> rdp;
  arr<float> rdt;
  float rn;
  float v1;
  float v2;
  float v3;
  float xj;
  float yj;

  quench_save() :
    amshu(fem::float0),
    bb(fem::float0),
    bbx(fem::float0),
    bby(fem::float0),
    de(fem::float0),
    dp(fem::float0),
    dp1(fem::float0),
    dp2(fem::float0),
    dp3(fem::float0),
    dphi(fem::float0),
    drr(fem::float0),
    dx(fem::float0),
    dy(fem::float0),
    ershu(fem::float0),
    i(fem::int0),
    i2(fem::int0),
    isg(fem::int0),
    j2(fem::int0),
    jp(fem::int0),
    jt(fem::int0),
    kp(fem::int0),
    kt(fem::int0),
    lq(fem::int0),
    lqp(dimension(300), fem::fill0),
    lqt(dimension(300), fem::fill0),
    mp(fem::int0),
    mt(fem::int0),
    nq(fem::int0),
    phi(fem::float0),
    phip(fem::float0),
    phiq(fem::float0),
    phit(fem::float0),
    prshu(fem::float0),
    ptjet0(fem::float0),
    ptot(fem::float0),
    r0(fem::float0),
    rd(fem::float0),
    rd0(fem::float0),
    rdp(dimension(300), fem::fill0),
    rdt(dimension(300), fem::fill0),
    rn(fem::float0),
    v1(fem::float0),
    v2(fem::float0),
    v3(fem::float0),
    xj(fem::float0),
    yj(fem::float0)
  {}
};

void
quench(
  common& cmn,
  int const& jpjt,
  int const& ntp)
{
  FEM_CMN_SVE(quench);
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_ref<int> npj(cmn.npj, dimension(300));
  arr_ref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_ref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_ref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_ref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_ref<float, 2> pjpe(cmn.pjpe, dimension(300, 500));
  arr_ref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_ref<int> ntj(cmn.ntj, dimension(300));
  arr_ref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_ref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_ref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_ref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_ref<float, 2> pjte(cmn.pjte, dimension(300, 500));
  arr_ref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  const int maxstr = 150001;
  arr_cref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_cref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_cref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_ref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_ref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_ref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_ref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_cref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  arr_cref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_cref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  int& nseed = cmn.nseed;
  //
  float& amshu = sve.amshu;
  float& bb = sve.bb;
  float& bbx = sve.bbx;
  float& bby = sve.bby;
  float& de = sve.de;
  float& dp = sve.dp;
  float& dp1 = sve.dp1;
  float& dp2 = sve.dp2;
  float& dp3 = sve.dp3;
  float& dphi = sve.dphi;
  float& drr = sve.drr;
  float& dx = sve.dx;
  float& dy = sve.dy;
  float& ershu = sve.ershu;
  int& i = sve.i;
  int& i2 = sve.i2;
  int& isg = sve.isg;
  int& j2 = sve.j2;
  int& jp = sve.jp;
  int& jt = sve.jt;
  int& kp = sve.kp;
  int& kt = sve.kt;
  int& lq = sve.lq;
  arr_ref<int> lqp(sve.lqp, dimension(300));
  arr_ref<int> lqt(sve.lqt, dimension(300));
  int& mp = sve.mp;
  int& mt = sve.mt;
  int& nq = sve.nq;
  float& phi = sve.phi;
  float& phip = sve.phip;
  float& phiq = sve.phiq;
  float& phit = sve.phit;
  float& prshu = sve.prshu;
  float& ptjet0 = sve.ptjet0;
  float& ptot = sve.ptot;
  float& r0 = sve.r0;
  float& rd = sve.rd;
  float& rd0 = sve.rd0;
  arr_ref<float> rdp(sve.rdp, dimension(300));
  arr_ref<float> rdt(sve.rdt, dimension(300));
  float& rn = sve.rn;
  float& v1 = sve.v1;
  float& v2 = sve.v2;
  float& v3 = sve.v3;
  float& xj = sve.xj;
  float& yj = sve.yj;
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HPARNT/
  //C
  //Cc      SAVE /HJJET1/
  //Cc      SAVE /HJJET2/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /RNDF77/
  //C
  //C     Uzhi:
  bb = hint1(19);
  phi = hint1(20);
  bbx = bb * fem::cos(phi);
  bby = bb * fem::sin(phi);
  //C
  if (ntp == 2) {
    goto statement_400;
  }
  if (ntp == 3) {
    goto statement_2000;
  }
  //C*******************************************************
  //C Jet interaction for proj jet in the direction PHIP
  //C******************************************************
  //C
  if (nfp(jpjt, 7) != 1) {
    return;
  }
  //C
  jp = jpjt;
  FEM_DO_SAFE(i, 1, npj(jp)) {
    ptjet0 = fem::sqrt(fem::pow2(pjpx(jp, i)) + fem::pow2(pjpy(jp, i)));
    if (ptjet0 <= hipr1(11)) {
      goto statement_290;
    }
    ptot = fem::sqrt(ptjet0 * ptjet0 + fem::pow2(pjpz(jp, i)));
    if (ptot < hipr1(8)) {
      goto statement_290;
    }
    phip = ulangl(pjpx(jp, i), pjpy(jp, i));
    //C******* find the wounded proj which can interact with jet***
    kp = 0;
    FEM_DO_SAFE(i2, 1, ihnt2(1)) {
      if (nfp(i2, 5) != 3 || i2 == jp) {
        goto statement_100;
      }
      dx = yp(1, i2) - yp(1, jp);
      dy = yp(2, i2) - yp(2, jp);
      phi = ulangl(dx, dy);
      dphi = fem::abs(phi - phip);
      //C     Uzhi:
      if (dphi >= hipr1(40)) {
        dphi = 2.f * hipr1(40) - dphi;
      }
      if (dphi >= hipr1(40) / 2.0f) {
        goto statement_100;
      }
      rd0 = fem::sqrt(dx * dx + dy * dy);
      if (rd0 * fem::sin(dphi) > hipr1(12)) {
        goto statement_100;
      }
      kp++;
      lqp(kp) = i2;
      rdp(kp) = fem::cos(dphi) * rd0;
      statement_100:;
    }
    //C*******        rearrange according decending rd************
    FEM_DO_SAFE(i2, 1, kp - 1) {
      FEM_DO_SAFE(j2, i2 + 1, kp) {
        if (rdp(i2) < rdp(j2)) {
          goto statement_110;
        }
        rd = rdp(i2);
        lq = lqp(i2);
        rdp(i2) = rdp(j2);
        lqp(i2) = lqp(j2);
        rdp(j2) = rd;
        lqp(j2) = lq;
        statement_110:;
      }
    }
    //C****** find wounded targ which can interact with jet********
    kt = 0;
    FEM_DO_SAFE(i2, 1, ihnt2(3)) {
      if (nft(i2, 5) != 3) {
        goto statement_120;
      }
      dx = yt(1, i2) - yp(1, jp) - bbx;
      dy = yt(2, i2) - yp(2, jp) - bby;
      phi = ulangl(dx, dy);
      dphi = fem::abs(phi - phip);
      //C     Uzhi:
      if (dphi >= hipr1(40)) {
        dphi = 2.f * hipr1(40) - dphi;
      }
      if (dphi > hipr1(40) / 2.0f) {
        goto statement_120;
      }
      rd0 = fem::sqrt(dx * dx + dy * dy);
      if (rd0 * fem::sin(dphi) > hipr1(12)) {
        goto statement_120;
      }
      kt++;
      lqt(kt) = i2;
      rdt(kt) = fem::cos(dphi) * rd0;
      statement_120:;
    }
    //C*******        rearrange according decending rd************
    FEM_DO_SAFE(i2, 1, kt - 1) {
      FEM_DO_SAFE(j2, i2 + 1, kt) {
        if (rdt(i2) < rdt(j2)) {
          goto statement_130;
        }
        rd = rdt(i2);
        lq = lqt(i2);
        rdt(i2) = rdt(j2);
        lqt(i2) = lqt(j2);
        rdt(j2) = rd;
        lqt(j2) = lq;
        statement_130:;
      }
    }
    //C
    mp = 0;
    mt = 0;
    r0 = 0.0f;
    nq = 0;
    dp = 0.0f;
    ptot = fem::sqrt(fem::pow2(pjpx(jp, i)) + fem::pow2(pjpy(jp,
      i)) + fem::pow2(pjpz(jp, i)));
    v1 = pjpx(jp, i) / ptot;
    v2 = pjpy(jp, i) / ptot;
    v3 = pjpz(jp, i) / ptot;
    //C
    statement_200:
    rn = ranart(nseed);
    statement_210:
    if (mt >= kt && mp >= kp) {
      goto statement_290;
    }
    if (mt >= kt) {
      goto statement_220;
    }
    if (mp >= kp) {
      goto statement_240;
    }
    if (rdp(mp + 1) > rdt(mt + 1)) {
      goto statement_240;
    }
    statement_220:
    mp++;
    drr = rdp(mp) - r0;
    if (rn >= 1.0f - fem::exp(-drr / hipr1(13))) {
      goto statement_210;
    }
    dp = drr * hipr1(14);
    if (kfpj(jp, i) != 21) {
      dp = 0.5f * dp;
    }
    //C        ********string tension of quark jet is 0.5 of gluon's
    if (dp <= 0.2f) {
      goto statement_210;
    }
    if (ptot <= 0.4f) {
      goto statement_290;
    }
    if (ptot <= dp) {
      dp = ptot - 0.2f;
    }
    de = dp;
    //C
    if (kfpj(jp, i) != 21) {
      prshu = fem::pow2(pp(lqp(mp), 1)) + fem::pow2(pp(lqp(mp), 2)) +
        fem::pow2(pp(lqp(mp), 3));
      de = fem::sqrt(fem::pow2(pjpm(jp, i)) + fem::pow2(ptot)) -
        fem::sqrt(fem::pow2(pjpm(jp, i)) + fem::pow2((ptot - dp)));
      ershu = fem::pow2((pp(lqp(mp), 4) + de - dp));
      amshu = ershu - prshu;
      if (amshu < hipr1(1) * hipr1(1)) {
        goto statement_210;
      }
      pp(lqp(mp), 4) = fem::sqrt(ershu);
      pp(lqp(mp), 5) = fem::sqrt(amshu);
    }
    //C                ********reshuffle the energy when jet has mass
    r0 = rdp(mp);
    dp1 = dp * v1;
    dp2 = dp * v2;
    dp3 = dp * v3;
    //C                ********momentum and energy transfer from jet
    //C
    npj(lqp(mp))++;
    kfpj(lqp(mp), npj(lqp(mp))) = 21;
    pjpx(lqp(mp), npj(lqp(mp))) = dp1;
    pjpy(lqp(mp), npj(lqp(mp))) = dp2;
    pjpz(lqp(mp), npj(lqp(mp))) = dp3;
    pjpe(lqp(mp), npj(lqp(mp))) = dp;
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0f;
    goto statement_260;
    //C
    statement_240:
    mt++;
    drr = rdt(mt) - r0;
    if (rn >= 1.0f - fem::exp(-drr / hipr1(13))) {
      goto statement_210;
    }
    dp = drr * hipr1(14);
    if (dp <= 0.2f) {
      goto statement_210;
    }
    if (ptot <= 0.4f) {
      goto statement_290;
    }
    if (ptot <= dp) {
      dp = ptot - 0.2f;
    }
    de = dp;
    //C
    if (kfpj(jp, i) != 21) {
      prshu = fem::pow2(pt(lqt(mt), 1)) + fem::pow2(pt(lqt(mt), 2)) +
        fem::pow2(pt(lqt(mt), 3));
      de = fem::sqrt(fem::pow2(pjpm(jp, i)) + fem::pow2(ptot)) -
        fem::sqrt(fem::pow2(pjpm(jp, i)) + fem::pow2((ptot - dp)));
      ershu = fem::pow2((pt(lqt(mt), 4) + de - dp));
      amshu = ershu - prshu;
      if (amshu < hipr1(1) * hipr1(1)) {
        goto statement_210;
      }
      pt(lqt(mt), 4) = fem::sqrt(ershu);
      pt(lqt(mt), 5) = fem::sqrt(amshu);
    }
    //C                ********reshuffle the energy when jet has mass
    //C
    r0 = rdt(mt);
    dp1 = dp * v1;
    dp2 = dp * v2;
    dp3 = dp * v3;
    //C                ********momentum and energy transfer from jet
    ntj(lqt(mt))++;
    kftj(lqt(mt), ntj(lqt(mt))) = 21;
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1;
    pjty(lqt(mt), ntj(lqt(mt))) = dp2;
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3;
    pjte(lqt(mt), ntj(lqt(mt))) = dp;
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0f;
    //C
    statement_260:
    pjpx(jp, i) = (ptot - dp) * v1;
    pjpy(jp, i) = (ptot - dp) * v2;
    pjpz(jp, i) = (ptot - dp) * v3;
    pjpe(jp, i) = pjpe(jp, i) - de;
    //C
    ptot = ptot - dp;
    nq++;
    goto statement_200;
    statement_290:;
  }
  //C
  return;
  //C
  //C*******************************************************
  //C Jet interaction for target jet in the direction PHIT
  //C******************************************************
  //C
  //C******* find the wounded proj which can interact with jet***
  //C
  statement_400:
  if (nft(jpjt, 7) != 1) {
    return;
  }
  jt = jpjt;
  FEM_DO_SAFE(i, 1, ntj(jt)) {
    ptjet0 = fem::sqrt(fem::pow2(pjtx(jt, i)) + fem::pow2(pjty(jt, i)));
    if (ptjet0 <= hipr1(11)) {
      goto statement_690;
    }
    ptot = fem::sqrt(ptjet0 * ptjet0 + fem::pow2(pjtz(jt, i)));
    if (ptot < hipr1(8)) {
      goto statement_690;
    }
    phit = ulangl(pjtx(jt, i), pjty(jt, i));
    kp = 0;
    FEM_DO_SAFE(i2, 1, ihnt2(1)) {
      if (nfp(i2, 5) != 3) {
        goto statement_500;
      }
      dx = yp(1, i2) + bbx - yt(1, jt);
      dy = yp(2, i2) + bby - yt(2, jt);
      phi = ulangl(dx, dy);
      dphi = fem::abs(phi - phit);
      //C     Uzhi:
      if (dphi >= hipr1(40)) {
        dphi = 2.f * hipr1(40) - dphi;
      }
      if (dphi > hipr1(40) / 2.0f) {
        goto statement_500;
      }
      rd0 = fem::sqrt(dx * dx + dy * dy);
      if (rd0 * fem::sin(dphi) > hipr1(12)) {
        goto statement_500;
      }
      kp++;
      lqp(kp) = i2;
      rdp(kp) = fem::cos(dphi) * rd0;
      statement_500:;
    }
    //C*******        rearrange according to decending rd************
    FEM_DO_SAFE(i2, 1, kp - 1) {
      FEM_DO_SAFE(j2, i2 + 1, kp) {
        if (rdp(i2) < rdp(j2)) {
          goto statement_510;
        }
        rd = rdp(i2);
        lq = lqp(i2);
        rdp(i2) = rdp(j2);
        lqp(i2) = lqp(j2);
        rdp(j2) = rd;
        lqp(j2) = lq;
        statement_510:;
      }
    }
    //C****** find wounded targ which can interact with jet********
    kt = 0;
    FEM_DO_SAFE(i2, 1, ihnt2(3)) {
      if (nft(i2, 5) != 3 || i2 == jt) {
        goto statement_520;
      }
      dx = yt(1, i2) - yt(1, jt);
      dy = yt(2, i2) - yt(2, jt);
      phi = ulangl(dx, dy);
      dphi = fem::abs(phi - phit);
      //C     Uzhi:
      if (dphi >= hipr1(40)) {
        dphi = 2.f * hipr1(40) - dphi;
      }
      if (dphi > hipr1(40) / 2.0f) {
        goto statement_520;
      }
      rd0 = fem::sqrt(dx * dx + dy * dy);
      if (rd0 * fem::sin(dphi) > hipr1(12)) {
        goto statement_520;
      }
      kt++;
      lqt(kt) = i2;
      rdt(kt) = fem::cos(dphi) * rd0;
      statement_520:;
    }
    //C*******        rearrange according to decending rd************
    FEM_DO_SAFE(i2, 1, kt - 1) {
      FEM_DO_SAFE(j2, i2 + 1, kt) {
        if (rdt(i2) < rdt(j2)) {
          goto statement_530;
        }
        rd = rdt(i2);
        lq = lqt(i2);
        rdt(i2) = rdt(j2);
        lqt(i2) = lqt(j2);
        rdt(j2) = rd;
        lqt(j2) = lq;
        statement_530:;
      }
    }
    //C
    mp = 0;
    mt = 0;
    nq = 0;
    dp = 0.0f;
    r0 = 0.0f;
    ptot = fem::sqrt(fem::pow2(pjtx(jt, i)) + fem::pow2(pjty(jt,
      i)) + fem::pow2(pjtz(jt, i)));
    v1 = pjtx(jt, i) / ptot;
    v2 = pjty(jt, i) / ptot;
    v3 = pjtz(jt, i) / ptot;
    //C
    statement_600:
    rn = ranart(nseed);
    statement_610:
    if (mt >= kt && mp >= kp) {
      goto statement_690;
    }
    if (mt >= kt) {
      goto statement_620;
    }
    if (mp >= kp) {
      goto statement_640;
    }
    if (rdp(mp + 1) > rdt(mt + 1)) {
      goto statement_640;
    }
    statement_620:
    mp++;
    drr = rdp(mp) - r0;
    if (rn >= 1.0f - fem::exp(-drr / hipr1(13))) {
      goto statement_610;
    }
    dp = drr * hipr1(14);
    if (kftj(jt, i) != 21) {
      dp = 0.5f * dp;
    }
    //C        ********string tension of quark jet is 0.5 of gluon's
    if (dp <= 0.2f) {
      goto statement_610;
    }
    if (ptot <= 0.4f) {
      goto statement_690;
    }
    if (ptot <= dp) {
      dp = ptot - 0.2f;
    }
    de = dp;
    //C
    if (kftj(jt, i) != 21) {
      prshu = fem::pow2(pp(lqp(mp), 1)) + fem::pow2(pp(lqp(mp), 2)) +
        fem::pow2(pp(lqp(mp), 3));
      de = fem::sqrt(fem::pow2(pjtm(jt, i)) + fem::pow2(ptot)) -
        fem::sqrt(fem::pow2(pjtm(jt, i)) + fem::pow2((ptot - dp)));
      ershu = fem::pow2((pp(lqp(mp), 4) + de - dp));
      amshu = ershu - prshu;
      if (amshu < hipr1(1) * hipr1(1)) {
        goto statement_610;
      }
      pp(lqp(mp), 4) = fem::sqrt(ershu);
      pp(lqp(mp), 5) = fem::sqrt(amshu);
    }
    //C                ********reshuffle the energy when jet has mass
    //C
    r0 = rdp(mp);
    dp1 = dp * v1;
    dp2 = dp * v2;
    dp3 = dp * v3;
    //C                ********momentum and energy transfer from jet
    npj(lqp(mp))++;
    kfpj(lqp(mp), npj(lqp(mp))) = 21;
    pjpx(lqp(mp), npj(lqp(mp))) = dp1;
    pjpy(lqp(mp), npj(lqp(mp))) = dp2;
    pjpz(lqp(mp), npj(lqp(mp))) = dp3;
    pjpe(lqp(mp), npj(lqp(mp))) = dp;
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0f;
    //C
    goto statement_660;
    //C
    statement_640:
    mt++;
    drr = rdt(mt) - r0;
    if (rn >= 1.0f - fem::exp(-drr / hipr1(13))) {
      goto statement_610;
    }
    dp = drr * hipr1(14);
    if (dp <= 0.2f) {
      goto statement_610;
    }
    if (ptot <= 0.4f) {
      goto statement_690;
    }
    if (ptot <= dp) {
      dp = ptot - 0.2f;
    }
    de = dp;
    //C
    if (kftj(jt, i) != 21) {
      prshu = fem::pow2(pt(lqt(mt), 1)) + fem::pow2(pt(lqt(mt), 2)) +
        fem::pow2(pt(lqt(mt), 3));
      de = fem::sqrt(fem::pow2(pjtm(jt, i)) + fem::pow2(ptot)) -
        fem::sqrt(fem::pow2(pjtm(jt, i)) + fem::pow2((ptot - dp)));
      ershu = fem::pow2((pt(lqt(mt), 4) + de - dp));
      amshu = ershu - prshu;
      if (amshu < hipr1(1) * hipr1(1)) {
        goto statement_610;
      }
      pt(lqt(mt), 4) = fem::sqrt(ershu);
      pt(lqt(mt), 5) = fem::sqrt(amshu);
    }
    //C                ********reshuffle the energy when jet has mass
    //C
    r0 = rdt(mt);
    dp1 = dp * v1;
    dp2 = dp * v2;
    dp3 = dp * v3;
    //C                ********momentum and energy transfer from jet
    ntj(lqt(mt))++;
    kftj(lqt(mt), ntj(lqt(mt))) = 21;
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1;
    pjty(lqt(mt), ntj(lqt(mt))) = dp2;
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3;
    pjte(lqt(mt), ntj(lqt(mt))) = dp;
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0f;
    //C
    statement_660:
    pjtx(jt, i) = (ptot - dp) * v1;
    pjty(jt, i) = (ptot - dp) * v2;
    pjtz(jt, i) = (ptot - dp) * v3;
    pjte(jt, i) = pjte(jt, i) - de;
    //C
    ptot = ptot - dp;
    nq++;
    goto statement_600;
    statement_690:;
  }
  return;
  //C********************************************************
  //C        Q-QBAR jet interaction
  //C********************************************************
  statement_2000:
  isg = jpjt;
  if (iasg(isg, 3) != 1) {
    return;
  }
  //C
  jp = iasg(isg, 1);
  jt = iasg(isg, 2);
  xj = (yp(1, jp) + bbx + yt(1, jt)) / 2.0f;
  yj = (yp(2, jp) + bby + yt(2, jt)) / 2.0f;
  FEM_DO_SAFE(i, 1, njsg(isg)) {
    ptjet0 = fem::sqrt(fem::pow2(pxsg(isg, i)) + fem::pow2(pysg(isg, i)));
    if (ptjet0 <= hipr1(11) || pesg(isg, i) < hipr1(1)) {
      goto statement_2690;
    }
    ptot = fem::sqrt(ptjet0 * ptjet0 + fem::pow2(pzsg(isg, i)));
    if (ptot < fem::max(hipr1(1), hipr1(8))) {
      goto statement_2690;
    }
    phiq = ulangl(pxsg(isg, i), pysg(isg, i));
    kp = 0;
    FEM_DO_SAFE(i2, 1, ihnt2(1)) {
      if (nfp(i2, 5) != 3 || i2 == jp) {
        goto statement_2500;
      }
      dx = yp(1, i2) + bbx - xj;
      dy = yp(2, i2) + bby - yj;
      phi = ulangl(dx, dy);
      dphi = fem::abs(phi - phiq);
      //C     Uzhi:
      if (dphi >= hipr1(40)) {
        dphi = 2.f * hipr1(40) - dphi;
      }
      if (dphi > hipr1(40) / 2.0f) {
        goto statement_2500;
      }
      rd0 = fem::sqrt(dx * dx + dy * dy);
      if (rd0 * fem::sin(dphi) > hipr1(12)) {
        goto statement_2500;
      }
      kp++;
      lqp(kp) = i2;
      rdp(kp) = fem::cos(dphi) * rd0;
      statement_2500:;
    }
    //C*******        rearrange according to decending rd************
    FEM_DO_SAFE(i2, 1, kp - 1) {
      FEM_DO_SAFE(j2, i2 + 1, kp) {
        if (rdp(i2) < rdp(j2)) {
          goto statement_2510;
        }
        rd = rdp(i2);
        lq = lqp(i2);
        rdp(i2) = rdp(j2);
        lqp(i2) = lqp(j2);
        rdp(j2) = rd;
        lqp(j2) = lq;
        statement_2510:;
      }
    }
    //C****** find wounded targ which can interact with jet********
    kt = 0;
    FEM_DO_SAFE(i2, 1, ihnt2(3)) {
      if (nft(i2, 5) != 3 || i2 == jt) {
        goto statement_2520;
      }
      dx = yt(1, i2) - xj;
      dy = yt(2, i2) - yj;
      phi = ulangl(dx, dy);
      dphi = fem::abs(phi - phiq);
      //C     Uzhi:
      if (dphi >= hipr1(40)) {
        dphi = 2.f * hipr1(40) - dphi;
      }
      if (dphi > hipr1(40) / 2.0f) {
        goto statement_2520;
      }
      rd0 = fem::sqrt(dx * dx + dy * dy);
      if (rd0 * fem::sin(dphi) > hipr1(12)) {
        goto statement_2520;
      }
      kt++;
      lqt(kt) = i2;
      rdt(kt) = fem::cos(dphi) * rd0;
      statement_2520:;
    }
    //C*******        rearrange according to decending rd************
    FEM_DO_SAFE(i2, 1, kt - 1) {
      FEM_DO_SAFE(j2, i2 + 1, kt) {
        if (rdt(i2) < rdt(j2)) {
          goto statement_2530;
        }
        rd = rdt(i2);
        lq = lqt(i2);
        rdt(i2) = rdt(j2);
        lqt(i2) = lqt(j2);
        rdt(j2) = rd;
        lqt(j2) = lq;
        statement_2530:;
      }
    }
    //C
    mp = 0;
    mt = 0;
    nq = 0;
    dp = 0.0f;
    r0 = 0.0f;
    ptot = fem::sqrt(fem::pow2(pxsg(isg, i)) + fem::pow2(pysg(isg,
      i)) + fem::pow2(pzsg(isg, i)));
    v1 = pxsg(isg, i) / ptot;
    v2 = pysg(isg, i) / ptot;
    v3 = pzsg(isg, i) / ptot;
    //C
    statement_2600:
    rn = ranart(nseed);
    statement_2610:
    if (mt >= kt && mp >= kp) {
      goto statement_2690;
    }
    if (mt >= kt) {
      goto statement_2620;
    }
    if (mp >= kp) {
      goto statement_2640;
    }
    if (rdp(mp + 1) > rdt(mt + 1)) {
      goto statement_2640;
    }
    statement_2620:
    mp++;
    drr = rdp(mp) - r0;
    if (rn >= 1.0f - fem::exp(-drr / hipr1(13))) {
      goto statement_2610;
    }
    dp = drr * hipr1(14) / 2.0f;
    if (dp <= 0.2f) {
      goto statement_2610;
    }
    if (ptot <= 0.4f) {
      goto statement_2690;
    }
    if (ptot <= dp) {
      dp = ptot - 0.2f;
    }
    de = dp;
    //C
    if (k2sg(isg, i) != 21) {
      if (ptot < dp + hipr1(1)) {
        goto statement_2690;
      }
      prshu = fem::pow2(pp(lqp(mp), 1)) + fem::pow2(pp(lqp(mp), 2)) +
        fem::pow2(pp(lqp(mp), 3));
      de = fem::sqrt(fem::pow2(pmsg(isg, i)) + fem::pow2(ptot)) -
        fem::sqrt(fem::pow2(pmsg(isg, i)) + fem::pow2((ptot - dp)));
      ershu = fem::pow2((pp(lqp(mp), 4) + de - dp));
      amshu = ershu - prshu;
      if (amshu < hipr1(1) * hipr1(1)) {
        goto statement_2610;
      }
      pp(lqp(mp), 4) = fem::sqrt(ershu);
      pp(lqp(mp), 5) = fem::sqrt(amshu);
    }
    //C                ********reshuffle the energy when jet has mass
    //C
    r0 = rdp(mp);
    dp1 = dp * v1;
    dp2 = dp * v2;
    dp3 = dp * v3;
    //C                ********momentum and energy transfer from jet
    npj(lqp(mp))++;
    kfpj(lqp(mp), npj(lqp(mp))) = 21;
    pjpx(lqp(mp), npj(lqp(mp))) = dp1;
    pjpy(lqp(mp), npj(lqp(mp))) = dp2;
    pjpz(lqp(mp), npj(lqp(mp))) = dp3;
    pjpe(lqp(mp), npj(lqp(mp))) = dp;
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0f;
    //C
    goto statement_2660;
    //C
    statement_2640:
    mt++;
    drr = rdt(mt) - r0;
    if (rn >= 1.0f - fem::exp(-drr / hipr1(13))) {
      goto statement_2610;
    }
    dp = drr * hipr1(14);
    if (dp <= 0.2f) {
      goto statement_2610;
    }
    if (ptot <= 0.4f) {
      goto statement_2690;
    }
    if (ptot <= dp) {
      dp = ptot - 0.2f;
    }
    de = dp;
    //C
    if (k2sg(isg, i) != 21) {
      if (ptot < dp + hipr1(1)) {
        goto statement_2690;
      }
      prshu = fem::pow2(pt(lqt(mt), 1)) + fem::pow2(pt(lqt(mt), 2)) +
        fem::pow2(pt(lqt(mt), 3));
      de = fem::sqrt(fem::pow2(pmsg(isg, i)) + fem::pow2(ptot)) -
        fem::sqrt(fem::pow2(pmsg(isg, i)) + fem::pow2((ptot - dp)));
      ershu = fem::pow2((pt(lqt(mt), 4) + de - dp));
      amshu = ershu - prshu;
      if (amshu < hipr1(1) * hipr1(1)) {
        goto statement_2610;
      }
      pt(lqt(mt), 4) = fem::sqrt(ershu);
      pt(lqt(mt), 5) = fem::sqrt(amshu);
    }
    //C               ********reshuffle the energy when jet has mass
    //C
    r0 = rdt(mt);
    dp1 = dp * v1;
    dp2 = dp * v2;
    dp3 = dp * v3;
    //C                ********momentum and energy transfer from jet
    ntj(lqt(mt))++;
    kftj(lqt(mt), ntj(lqt(mt))) = 21;
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1;
    pjty(lqt(mt), ntj(lqt(mt))) = dp2;
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3;
    pjte(lqt(mt), ntj(lqt(mt))) = dp;
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0f;
    //C
    statement_2660:
    pxsg(isg, i) = (ptot - dp) * v1;
    pysg(isg, i) = (ptot - dp) * v2;
    pzsg(isg, i) = (ptot - dp) * v3;
    pesg(isg, i) = pesg(isg, i) - de;
    //C
    ptot = ptot - dp;
    nq++;
    goto statement_2600;
    statement_2690:;
  }
}

struct ar3jet_save
{
  float a;
  float c;
  float d;
  float exp1;
  float exp3;
  float fg;
  int neg;
  int ntry;
  float sm1;
  float sm3;
  float x2;
  float xt2;
  float xt2m;
  float y;
  float yma;
  float ymax;

  ar3jet_save() :
    a(fem::float0),
    c(fem::float0),
    d(fem::float0),
    exp1(fem::float0),
    exp3(fem::float0),
    fg(fem::float0),
    neg(fem::int0),
    ntry(fem::int0),
    sm1(fem::float0),
    sm3(fem::float0),
    x2(fem::float0),
    xt2(fem::float0),
    xt2m(fem::float0),
    y(fem::float0),
    yma(fem::float0),
    ymax(fem::float0)
  {}
};

void
ar3jet(
  common& cmn,
  float const& s,
  float& x1,
  float& x3,
  int const& jl)
{
  FEM_CMN_SVE(ar3jet);
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_cref<float, 2> p(cmn.p, dimension(9000, 5));
  int& nseed = cmn.nseed;
  //
  float& a = sve.a;
  float& c = sve.c;
  float& d = sve.d;
  float& exp1 = sve.exp1;
  float& exp3 = sve.exp3;
  float& fg = sve.fg;
  int& neg = sve.neg;
  int& ntry = sve.ntry;
  float& sm1 = sve.sm1;
  float& sm3 = sve.sm3;
  float& x2 = sve.x2;
  float& xt2 = sve.xt2;
  float& xt2m = sve.xt2m;
  float& y = sve.y;
  float& yma = sve.yma;
  float& ymax = sve.ymax;
  //C
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /RNDF77/
  //C
  c = 1.f / 3.f;
  if (k(jl, 2) != 21 && k(jl + 1, 2) != 21) {
    c = 8.f / 27.f;
  }
  exp1 = 3;
  exp3 = 3;
  if (k(jl, 2) != 21) {
    exp1 = 2;
  }
  if (k(jl + 1, 2) != 21) {
    exp3 = 2;
  }
  a = fem::pow2(0.24f) / s;
  yma = fem::alog(.5f / fem::sqrt(a) + fem::sqrt(.25f / a - 1));
  d = 4.f * c * yma;
  sm1 = fem::pow2(p(jl, 5)) / s;
  sm3 = fem::pow2(p(jl + 1, 5)) / s;
  xt2m = (1.f - 2.f * fem::sqrt(sm1) + sm1 - sm3) * (1.f - 2.f *
    fem::sqrt(sm3) - sm1 + sm3);
  xt2m = fem::min(.25f, xt2m);
  ntry = 0;
  statement_1:
  if (ntry == 5000) {
    x1 = .5f * (2.f * fem::sqrt(sm1) + 1.f + sm1 - sm3);
    x3 = .5f * (2.f * fem::sqrt(sm3) + 1.f - sm1 + sm3);
    return;
  }
  ntry++;
  //C
  xt2 = a * fem::pow((xt2m / a), (fem::pow(ranart(nseed), (1.f / d))));
  //C
  ymax = fem::alog(.5f / fem::sqrt(xt2) + fem::sqrt(.25f / xt2 - 1.f));
  y = (2.f * ranart(nseed) - 1.f) * ymax;
  x1 = 1.f - fem::sqrt(xt2) * fem::exp(y);
  x3 = 1.f - fem::sqrt(xt2) * fem::exp(-y);
  x2 = 2.f - x1 - x3;
  neg = 0;
  if (k(jl, 2) != 21 || k(jl + 1, 2) != 21) {
    if ((1.f - x1) * (1.f - x2) * (1.f - x3) - x2 * sm1 * (1.f -
        x1) - x2 * sm3 * (1.f - x3) <= 0.f || x1 <= 2.f * fem::sqrt(
        sm1) - sm1 + sm3 || x3 <= 2.f * fem::sqrt(sm3) - sm3 + sm1) {
      neg = 1;
    }
    x1 += sm1 - sm3;
    x3 = x3 - sm1 + sm3;
  }
  if (neg == 1) {
    goto statement_1;
  }
  //C
  fg = 2.f * ymax * c * (fem::pow(x1, exp1) + fem::pow(x3, exp3)) / d;
  xt2m = xt2;
  if (fg < ranart(nseed)) {
    goto statement_1;
  }
  //C
}

struct arorie_save
{
  float bet;
  float cbet;
  float del;
  float e1;
  float e3;
  float p1;
  float p3;
  float psi;
  float pt1;
  float pt3;
  float pz1;
  float pz3;
  float w;
  float x2;

  arorie_save() :
    bet(fem::float0),
    cbet(fem::float0),
    del(fem::float0),
    e1(fem::float0),
    e3(fem::float0),
    p1(fem::float0),
    p3(fem::float0),
    psi(fem::float0),
    pt1(fem::float0),
    pt3(fem::float0),
    pz1(fem::float0),
    pz3(fem::float0),
    w(fem::float0),
    x2(fem::float0)
  {}
};

//C*************************************************************
//C
void
arorie(
  common& cmn,
  float const& s,
  float const& x1,
  float const& x3,
  int const& jl)
{
  FEM_CMN_SVE(arorie);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  // COMMON lujets
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  //
  // SAVE
  float& bet = sve.bet;
  float& cbet = sve.cbet;
  float& del = sve.del;
  float& e1 = sve.e1;
  float& e3 = sve.e3;
  float& p1 = sve.p1;
  float& p3 = sve.p3;
  float& psi = sve.psi;
  float& pt1 = sve.pt1;
  float& pt3 = sve.pt3;
  float& pz1 = sve.pz1;
  float& pz3 = sve.pz3;
  float& w = sve.w;
  float& x2 = sve.x2;
  //
  //C
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /RNDF77/
  //C
  w = fem::sqrt(s);
  x2 = 2.f - x1 - x3;
  e1 = .5f * x1 * w;
  e3 = .5f * x3 * w;
  p1 = fem::sqrt(fem::pow2(e1) - fem::pow2(p(jl, 5)));
  p3 = fem::sqrt(fem::pow2(e3) - fem::pow2(p(jl + 1, 5)));
  cbet = 1.f;
  if (p1 > 0.f && p3 > 0.f) {
    cbet = (fem::pow2(p(jl, 5)) + fem::pow2(p(jl + 1, 5)) + 2.f *
      e1 * e3 - s * (1.f - x2)) / (2.f * p1 * p3);
  }
  if (fem::abs(cbet) > 1.0f) {
    cbet = fem::max(-1.f, fem::min(1.f, cbet));
  }
  bet = fem::acos(cbet);
  //C
  //C.....MINIMIZE PT1-SQUARED PLUS PT3-SQUARED.....
  if (p1 >= p3) {
    psi = .5f * ulangl(fem::pow2(p1) + fem::pow2(p3) * fem::cos(2.f * bet),
      -fem::pow2(p3) * fem::sin(2.f * bet));
    pt1 = p1 * fem::sin(psi);
    pz1 = p1 * fem::cos(psi);
    pt3 = p3 * fem::sin(psi + bet);
    pz3 = p3 * fem::cos(psi + bet);
  }
  else if (p3 > p1) {
    psi = .5f * ulangl(fem::pow2(p3) + fem::pow2(p1) * fem::cos(2.f * bet),
      -fem::pow2(p1) * fem::sin(2.f * bet));
    pt1 = p1 * fem::sin(bet + psi);
    pz1 = -p1 * fem::cos(bet + psi);
    pt3 = p3 * fem::sin(psi);
    pz3 = -p3 * fem::cos(psi);
  }
  //C
  del = 2.0f * hipr1(40) * ranart(cmn.nseed);
  p(jl, 4) = e1;
  p(jl, 1) = pt1 * fem::sin(del);
  p(jl, 2) = -pt1 * fem::cos(del);
  p(jl, 3) = pz1;
  p(jl + 2, 4) = e3;
  p(jl + 2, 1) = pt3 * fem::sin(del);
  p(jl + 2, 2) = -pt3 * fem::cos(del);
  p(jl + 2, 3) = pz3;
  p(jl + 1, 4) = w - e1 - e3;
  p(jl + 1, 1) = -p(jl, 1) - p(jl + 2, 1);
  p(jl + 1, 2) = -p(jl, 2) - p(jl + 2, 2);
  p(jl + 1, 3) = -p(jl, 3) - p(jl + 2, 3);
}

struct atrobo_save
{
  double dbep;
  double dbex;
  double dbey;
  double dbez;
  double dga;
  double dga2;
  double dgabep;
  arr<double> dp;
  int i;
  int j;
  arr<float> pv;
  arr<float, 2> rot;

  atrobo_save() :
    dbep(fem::double0),
    dbex(fem::double0),
    dbey(fem::double0),
    dbez(fem::double0),
    dga(fem::double0),
    dga2(fem::double0),
    dgabep(fem::double0),
    dp(dimension(4), fem::fill0),
    i(fem::int0),
    j(fem::int0),
    pv(dimension(3), fem::fill0),
    rot(dimension(3, 3), fem::fill0)
  {}
};

//C
//C*******************************************************************
//C        make  boost and rotation to entries from IMIN to IMAX
//C*******************************************************************
void
atrobo(
  common& cmn,
  float const& the,
  float const& phi,
  float const& bex,
  float const& bey,
  float const& bez,
  int const& imin,
  int const& imax,
  int& ierror)
{
  FEM_CMN_SVE(atrobo);
  // COMMON lujets
  int& n = static_cast<common_lujets&>(cmn).n;
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  //
  // SAVE
  double& dbep = sve.dbep;
  double& dbex = sve.dbex;
  double& dbey = sve.dbey;
  double& dbez = sve.dbez;
  double& dga = sve.dga;
  double& dga2 = sve.dga2;
  double& dgabep = sve.dgabep;
  arr_ref<double> dp(sve.dp, dimension(4));
  int& i = sve.i;
  int& j = sve.j;
  arr_ref<float> pv(sve.pv, dimension(3));
  arr_ref<float, 2> rot(sve.rot, dimension(3, 3));
  //
  //Cc      SAVE /LUJETS/
  ierror = 0;
  //C
  if (imin <= 0 || imax > n || imin > imax) {
    return;
  }
  //C
  if (fem::pow2(the) + fem::pow2(phi) > 1e-20f) {
    //C...ROTATE (TYPICALLY FROM Z AXIS TO DIRECTION THETA,PHI)
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
      //C**************           IF(MOD(K(I,1)/10000,10).GE.6) GOTO 120
      FEM_DO_SAFE(j, 1, 3) {
        pv(j) = p(i, j);
      }
      FEM_DO_SAFE(j, 1, 3) {
        p(i, j) = rot(j, 1) * pv(1) + rot(j, 2) * pv(2) + rot(j, 3) * pv(3);
      }
    }
  }
  //C
  if (fem::pow2(bex) + fem::pow2(bey) + fem::pow2(bez) > 1e-20f) {
    //C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
    dbex = fem::dble(bex);
    dbey = fem::dble(bey);
    dbez = fem::dble(bez);
    dga2 = 1e0 - fem::pow2(dbex) - fem::pow2(dbey) - fem::pow2(dbez);
    if (dga2 <= 0e0) {
      ierror = 1;
      return;
    }
    dga = 1e0 / fem::dsqrt(dga2);
    FEM_DO_SAFE(i, imin, imax) {
      //C*************           IF(MOD(K(I,1)/10000,10).GE.6) GOTO 140
      FEM_DO_SAFE(j, 1, 4) {
        dp(j) = fem::dble(p(i, j));
      }
      dbep = dbex * dp(1) + dbey * dp(2) + dbez * dp(3);
      dgabep = dga * (dga * dbep / (1e0 + dga) + dp(4));
      p(i, 1) = fem::sngl(dp(1) + dgabep * dbex);
      p(i, 2) = fem::sngl(dp(2) + dgabep * dbey);
      p(i, 3) = fem::sngl(dp(3) + dgabep * dbez);
      p(i, 4) = fem::sngl(dga * (dp(4) + dbep));
    }
  }
  //C
}

struct attrad_save
{
  float bex;
  float bey;
  float bez;
  float btt;
  float cth;
  float fmfact;
  int i;
  int imax;
  int imin;
  int j;
  int jl;
  int m;
  float p1;
  float p2;
  float p3;
  float p4;
  float pbt1;
  float pbt2;
  float pbt3;
  float pbt4;
  float phi;
  float ptg;
  float ptg1;
  float ptg2;
  float ptg3;
  float s;
  float sm;
  float theta;
  float wp;
  float x1;
  float x3;

  attrad_save() :
    bex(fem::float0),
    bey(fem::float0),
    bez(fem::float0),
    btt(fem::float0),
    cth(fem::float0),
    fmfact(fem::float0),
    i(fem::int0),
    imax(fem::int0),
    imin(fem::int0),
    j(fem::int0),
    jl(fem::int0),
    m(fem::int0),
    p1(fem::float0),
    p2(fem::float0),
    p3(fem::float0),
    p4(fem::float0),
    pbt1(fem::float0),
    pbt2(fem::float0),
    pbt3(fem::float0),
    pbt4(fem::float0),
    phi(fem::float0),
    ptg(fem::float0),
    ptg1(fem::float0),
    ptg2(fem::float0),
    ptg3(fem::float0),
    s(fem::float0),
    sm(fem::float0),
    theta(fem::float0),
    wp(fem::float0),
    x1(fem::float0),
    x3(fem::float0)
  {}
};

//C
//C****************************************************************
//C        conduct soft radiation according to dipole approxiamtion
//C****************************************************************
void
attrad(
  common& cmn,
  int& ierror)
{
  FEM_CMN_SVE(attrad);
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hidat(cmn.hidat, dimension(10));
  int& n = static_cast<common_lujets&>(cmn).n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  //
  float& bex = sve.bex;
  float& bey = sve.bey;
  float& bez = sve.bez;
  float& btt = sve.btt;
  float& cth = sve.cth;
  float& fmfact = sve.fmfact;
  int& i = sve.i;
  int& imax = sve.imax;
  int& imin = sve.imin;
  int& j = sve.j;
  int& jl = sve.jl;
  int& m = sve.m;
  float& p1 = sve.p1;
  float& p2 = sve.p2;
  float& p3 = sve.p3;
  float& p4 = sve.p4;
  float& pbt1 = sve.pbt1;
  float& pbt2 = sve.pbt2;
  float& pbt3 = sve.pbt3;
  float& pbt4 = sve.pbt4;
  float& phi = sve.phi;
  float& ptg = sve.ptg;
  float& ptg1 = sve.ptg1;
  float& ptg2 = sve.ptg2;
  float& ptg3 = sve.ptg3;
  float& s = sve.s;
  float& sm = sve.sm;
  float& theta = sve.theta;
  float& wp = sve.wp;
  float& x1 = sve.x1;
  float& x3 = sve.x3;
  //C
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HIJDAT/
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /RNDF77/
  ierror = 0;
  //C
  //C.....S INVARIANT MASS-SQUARED BETWEEN PARTONS I AND I+1......
  //C.....SM IS THE LARGEST MASS-SQUARED....
  //C
  statement_40:
  sm = 0.f;
  jl = 1;
  FEM_DO_SAFE(i, 1, n - 1) {
    s = 2.f * (p(i, 4) * p(i + 1, 4) - p(i, 1) * p(i + 1, 1) - p(i,
      2) * p(i + 1, 2) - p(i, 3) * p(i + 1, 3)) + fem::pow2(p(i,
      5)) + fem::pow2(p(i + 1, 5));
    if (s < 0.f) {
      s = 0.f;
    }
    wp = fem::sqrt(s) - 1.5f * (p(i, 5) + p(i + 1, 5));
    if (wp > sm) {
      pbt1 = p(i, 1) + p(i + 1, 1);
      pbt2 = p(i, 2) + p(i + 1, 2);
      pbt3 = p(i, 3) + p(i + 1, 3);
      pbt4 = p(i, 4) + p(i + 1, 4);
      btt = (fem::pow2(pbt1) + fem::pow2(pbt2) + fem::pow2(pbt3)) /
        fem::pow2(pbt4);
      if (btt >= 1.0f - 1.0e-10f) {
        goto statement_30;
      }
      if ((i != 1 || i != n - 1) && (k(i, 2) != 21 && k(i + 1, 2) != 21)) {
        goto statement_30;
      }
      jl = i;
      sm = wp;
    }
    statement_30:;
  }
  s = fem::pow2((sm + 1.5f * (p(jl, 5) + p(jl + 1, 5))));
  if (sm < hipr1(5)) {
    goto statement_2;
  }
  //C
  //C.....MAKE PLACE FOR ONE GLUON.....
  if (jl + 1 == n) {
    goto statement_190;
  }
  FEM_DOSTEP(j, n, jl + 2, -1) {
    k(j + 1, 1) = k(j, 1);
    k(j + 1, 2) = k(j, 2);
    FEM_DO_SAFE(m, 1, 5) {
      p(j + 1, m) = p(j, m);
    }
  }
  statement_190:
  n++;
  //C
  //C.....BOOST TO REST SYSTEM FOR PARTICLES JL AND JL+1.....
  p1 = p(jl, 1) + p(jl + 1, 1);
  p2 = p(jl, 2) + p(jl + 1, 2);
  p3 = p(jl, 3) + p(jl + 1, 3);
  p4 = p(jl, 4) + p(jl + 1, 4);
  bex = -p1 / p4;
  bey = -p2 / p4;
  bez = -p3 / p4;
  imin = jl;
  imax = jl + 1;
  atrobo(cmn, 0.f, 0.f, bex, bey, bez, imin, imax, ierror);
  if (ierror != 0) {
    return;
  }
  //C.....ROTATE TO Z-AXIS....
  cth = p(jl, 3) / fem::sqrt(fem::pow2(p(jl, 4)) - fem::pow2(p(jl, 5)));
  if (fem::abs(cth) > 1.0f) {
    cth = fem::max(-1.f, fem::min(1.f, cth));
  }
  theta = fem::acos(cth);
  phi = ulangl(p(jl, 1), p(jl, 2));
  atrobo(cmn, 0.f, -phi, 0.f, 0.f, 0.f, imin, imax, ierror);
  atrobo(cmn, -theta, 0.f, 0.f, 0.f, 0.f, imin, imax, ierror);
  //C
  //C.....CREATE ONE GLUON AND ORIENTATE.....
  //C
  statement_1:
  ar3jet(cmn, s, x1, x3, jl);
  arorie(cmn, s, x1, x3, jl);
  if (hidat(2) > 0.0f) {
    ptg1 = fem::sqrt(fem::pow2(p(jl, 1)) + fem::pow2(p(jl, 2)));
    ptg2 = fem::sqrt(fem::pow2(p(jl + 1, 1)) + fem::pow2(p(jl + 1, 2)));
    ptg3 = fem::sqrt(fem::pow2(p(jl + 2, 1)) + fem::pow2(p(jl + 2, 2)));
    ptg = fem::max(ptg1, ptg2, ptg3);
    if (ptg > hidat(2)) {
      fmfact = fem::exp(-(fem::pow2(ptg) - fem::pow2(hidat(2))) /
        fem::pow2(hipr1(2)));
      if (ranart(cmn.nseed) > fmfact) {
        goto statement_1;
      }
    }
  }
  //C.....ROTATE AND BOOST BACK.....
  imin = jl;
  imax = jl + 2;
  atrobo(cmn, theta, phi, -bex, -bey, -bez, imin, imax, ierror);
  if (ierror != 0) {
    return;
  }
  //C.....ENUMERATE THE GLUONS.....
  k(jl + 2, 1) = k(jl + 1, 1);
  k(jl + 2, 2) = k(jl + 1, 2);
  k(jl + 2, 3) = k(jl + 1, 3);
  k(jl + 2, 4) = k(jl + 1, 4);
  k(jl + 2, 5) = k(jl + 1, 5);
  p(jl + 2, 5) = p(jl + 1, 5);
  k(jl + 1, 1) = 2;
  k(jl + 1, 2) = 21;
  k(jl + 1, 3) = 0;
  k(jl + 1, 4) = 0;
  k(jl + 1, 5) = 0;
  p(jl + 1, 5) = 0.f;
  //C----THETA FUNCTION DAMPING OF THE EMITTED GLUONS. FOR HADRON-HADRON.
  //C----R0=VFR(2)
  //C              IF(VFR(2).GT.0.) THEN
  //C              PTG=SQRT(P(JL+1,1)**2+P(JL+1,2)**2)
  //C              PTGMAX=WSTRI/2.
  //C              DOPT=SQRT((4.*PAR(71)*VFR(2))/WSTRI)
  //C              PTOPT=(DOPT*WSTRI)/(2.*VFR(2))
  //C              IF(PTG.GT.PTOPT) IORDER=IORDER-1
  //C              IF(PTG.GT.PTOPT) GOTO 1
  //C              ENDIF
  //C-----
  if (sm >= hipr1(5)) {
    goto statement_40;
  }
  //C
  statement_2:
  k(1, 1) = 2;
  k(1, 3) = 0;
  k(1, 4) = 0;
  k(1, 5) = 0;
  k(n, 1) = 1;
  k(n, 3) = 0;
  k(n, 4) = 0;
  k(n, 5) = 0;
  //C
}

struct hijfrg_save
{
  float am1;
  float am2;
  float amt;
  float amt1;
  float amt2;
  float btz;
  float hdat20;
  float hpr150;
  int i;
  int i0;
  int iex;
  int ii;
  int isg;
  int j;
  int jetot;
  int jj;
  int kf1;
  int kf2;
  int kk1;
  float pb1;
  float pb2;
  float pb3;
  float pecm;
  float pmax1;
  float pmax2;
  float pmax3;
  float pq11;
  float pq12;
  float pq21;
  float pq22;
  float pzcm;

  hijfrg_save() :
    am1(fem::float0),
    am2(fem::float0),
    amt(fem::float0),
    amt1(fem::float0),
    amt2(fem::float0),
    btz(fem::float0),
    hdat20(fem::float0),
    hpr150(fem::float0),
    i(fem::int0),
    i0(fem::int0),
    iex(fem::int0),
    ii(fem::int0),
    isg(fem::int0),
    j(fem::int0),
    jetot(fem::int0),
    jj(fem::int0),
    kf1(fem::int0),
    kf2(fem::int0),
    kk1(fem::int0),
    pb1(fem::float0),
    pb2(fem::float0),
    pb3(fem::float0),
    pecm(fem::float0),
    pmax1(fem::float0),
    pmax2(fem::float0),
    pmax3(fem::float0),
    pq11(fem::float0),
    pq12(fem::float0),
    pq21(fem::float0),
    pq22(fem::float0),
    pzcm(fem::float0)
  {}
};

void
hijfrg(
  common& cmn,
  int const& jtp,
  int const& ntp,
  int& ierror)
{
  FEM_CMN_SVE(hijfrg);
  common_write write(cmn);
  arr_ref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_ref<float> hidat(cmn.hidat, dimension(10));
  arr_cref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_cref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  arr_cref<int> npj(cmn.npj, dimension(300));
  arr_cref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_cref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_cref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_cref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_cref<float, 2> pjpe(cmn.pjpe, dimension(300, 500));
  arr_cref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_cref<int> ntj(cmn.ntj, dimension(300));
  arr_cref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_cref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_cref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_cref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_cref<float, 2> pjte(cmn.pjte, dimension(300, 500));
  arr_cref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  const int maxstr = 150001;
  arr_cref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_cref<int, 2> k1sg(cmn.k1sg, dimension(maxstr, 100));
  arr_cref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_cref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_cref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_cref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_cref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_cref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  int& n = static_cast<common_lujets&>(cmn).n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<float, 2> v(cmn.v, dimension(9000, 5));
  int& nseed = cmn.nseed;
  int& isoft = cmn.isoft;
  int& isflag = cmn.isflag;
  //
  float& am1 = sve.am1;
  float& am2 = sve.am2;
  float& amt = sve.amt;
  float& amt1 = sve.amt1;
  float& amt2 = sve.amt2;
  float& btz = sve.btz;
  float& hdat20 = sve.hdat20;
  float& hpr150 = sve.hpr150;
  int& i = sve.i;
  int& i0 = sve.i0;
  int& iex = sve.iex;
  int& ii = sve.ii;
  int& isg = sve.isg;
  int& j = sve.j;
  int& jetot = sve.jetot;
  int& jj = sve.jj;
  int& kf1 = sve.kf1;
  int& kf2 = sve.kf2;
  int& kk1 = sve.kk1;
  float& pb1 = sve.pb1;
  float& pb2 = sve.pb2;
  float& pb3 = sve.pb3;
  float& pecm = sve.pecm;
  float& pmax1 = sve.pmax1;
  float& pmax2 = sve.pmax2;
  float& pmax3 = sve.pmax3;
  float& pq11 = sve.pq11;
  float& pq12 = sve.pq12;
  float& pq21 = sve.pq21;
  float& pq22 = sve.pq22;
  float& pzcm = sve.pzcm;
  //C        NTP=1, fragment proj string, NTP=2, targ string,
  //C       NTP=3, independent
  //C        strings from jets.  JTP is the line number of the string
  //C*******Fragment all leadng strings of proj and targ**************
  //C        IHNT2(1)=atomic #, IHNT2(2)=proton #(=-1 if anti-proton)  *
  //C******************************************************************
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HIJDAT/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /HJJET1/
  //Cc      SAVE /HJJET2/
  //C
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /LUDAT1/
  //Cc      SAVE /RNDF77/
  //Clin-4/11/01 soft:
  //Cc      SAVE /anim/
  //C
  //Cbz3/12/99
  //C.....set up fragmentation function according to the number of collisions
  //C.....a wounded nucleon has suffered
  //C        IF (NTP .EQ. 1) THEN
  //C           NCOLL = NFP(JTP, 11)
  //C        ELSE IF (NTP .EQ. 2) THEN
  //C           NCOLL = NFT(JTP, 11)
  //C        ELSE IF (NTP .EQ. 3) THEN
  //C           NCOLL = (NFP(IASG(JTP,1), 11) + NFT(IASG(JTP,2), 11)) / 2
  //C        END IF
  //C        IF (NCOLL .LE. 1) THEN
  //C           PARJ(5) = 0.5
  //C        ELSE IF (NCOLL .EQ. 2) THEN
  //C           PARJ(5) = 0.75
  //C        ELSE IF (NCOLL .EQ. 3) THEN
  //C           PARJ(5) = 1.17
  //C        ELSE IF (NCOLL .EQ. 4) THEN
  //C           PARJ(5) = 2.0
  //C        ELSE IF (NCOLL .EQ. 5) THEN
  //C           PARJ(5) = 4.5
  //C        ELSE IF (NCOLL .GE. 6) THEN
  //C           PARJ(5) = 49.5
  //C        END IF
  //C        PARJ(5) = 0.5
  //Cbz3/12/99 end
  //C
  ierror = 0;
  luedit(0);
  n = 0;
  //C                        ********initialize the document lines
  if (ntp == 3) {
    isg = jtp;
    n = njsg(isg);
    FEM_DO_SAFE(i, 1, njsg(isg)) {
      k(i, 1) = k1sg(isg, i);
      k(i, 2) = k2sg(isg, i);
      p(i, 1) = pxsg(isg, i);
      p(i, 2) = pysg(isg, i);
      p(i, 3) = pzsg(isg, i);
      p(i, 4) = pesg(isg, i);
      p(i, 5) = pmsg(isg, i);
    }
    //C
    //C                IF(IHPR2(1).GT.0) CALL ATTRAD(IERROR)
    //C                IF(IERROR.NE.0) RETURN
    //C                CALL LULIST(1)
    if (isoft != 2 || isflag != 0) {
      luexec();
    }
    return;
  }
  //C
  if (ntp == 2) {
    goto statement_200;
  }
  if (jtp > ihnt2(1)) {
    return;
  }
  if (nfp(jtp, 5) != 3 && nfp(jtp, 3) != 0 && npj(jtp) == 0 && nfp(jtp,
      10) == 0) {
    goto statement_1000;
  }
  if (nfp(jtp, 15) ==  - 1) {
    kf1 = nfp(jtp, 2);
    kf2 = nfp(jtp, 1);
    pq21 = pp(jtp, 6);
    pq22 = pp(jtp, 7);
    pq11 = pp(jtp, 8);
    pq12 = pp(jtp, 9);
    am1 = pp(jtp, 15);
    am2 = pp(jtp, 14);
  }
  else {
    kf1 = nfp(jtp, 1);
    kf2 = nfp(jtp, 2);
    pq21 = pp(jtp, 8);
    pq22 = pp(jtp, 9);
    pq11 = pp(jtp, 6);
    pq12 = pp(jtp, 7);
    am1 = pp(jtp, 14);
    am2 = pp(jtp, 15);
  }
  //C
  //C        ********for NFP(JTP,15)=-1 NFP(JTP,1) IS IN -Z DIRECTION
  pb1 = pq11 + pq21;
  pb2 = pq12 + pq22;
  pb3 = pp(jtp, 3);
  pecm = pp(jtp, 5);
  btz = pb3 / pp(jtp, 4);
  if ((fem::abs(pb1 - pp(jtp, 1)) > 0.01f || fem::abs(pb2 - pp(jtp,
      2)) > 0.01f) && ihpr2(10) != 0) {
    write(6, star), "  Pt of Q and QQ do not sum to the total", jtp,
      ntp, pq11, pq21, pb1, "*", pq12, pq22, pb2, "*", pp(jtp, 1), pp(jtp,
      2);
  }
  goto statement_300;
  //C
  statement_200:
  if (jtp > ihnt2(3)) {
    return;
  }
  if (nft(jtp, 5) != 3 && nft(jtp, 3) != 0 && ntj(jtp) == 0 && nft(jtp,
      10) == 0) {
    goto statement_1200;
  }
  if (nft(jtp, 15) == 1) {
    kf1 = nft(jtp, 1);
    kf2 = nft(jtp, 2);
    pq11 = pt(jtp, 6);
    pq12 = pt(jtp, 7);
    pq21 = pt(jtp, 8);
    pq22 = pt(jtp, 9);
    am1 = pt(jtp, 14);
    am2 = pt(jtp, 15);
  }
  else {
    kf1 = nft(jtp, 2);
    kf2 = nft(jtp, 1);
    pq11 = pt(jtp, 8);
    pq12 = pt(jtp, 9);
    pq21 = pt(jtp, 6);
    pq22 = pt(jtp, 7);
    am1 = pt(jtp, 15);
    am2 = pt(jtp, 14);
  }
  //C        ********for NFT(JTP,15)=1 NFT(JTP,1) IS IN +Z DIRECTION
  pb1 = pq11 + pq21;
  pb2 = pq12 + pq22;
  pb3 = pt(jtp, 3);
  pecm = pt(jtp, 5);
  btz = pb3 / pt(jtp, 4);
  //C
  if ((fem::abs(pb1 - pt(jtp, 1)) > 0.01f || fem::abs(pb2 - pt(jtp,
      2)) > 0.01f) && ihpr2(10) != 0) {
    write(6, star), "  Pt of Q and QQ do not sum to the total", jtp,
      ntp, pq11, pq21, pb1, "*", pq12, pq22, pb2, "*", pt(jtp, 1), pt(jtp,
      2);
  }
  statement_300:
  if (pecm < hipr1(1)) {
    ierror = 1;
    if (ihpr2(10) == 0) {
      return;
    }
    write(6, star), " ECM=", pecm, " energy of the string is too small";
    //Clin:
    write(6, star), "JTP,NTP,pq=", jtp, ntp, pq11, pq12, pq21, pq22;
    return;
  }
  amt = fem::pow2(pecm) + fem::pow2(pb1) + fem::pow2(pb2);
  amt1 = fem::pow2(am1) + fem::pow2(pq11) + fem::pow2(pq12);
  amt2 = fem::pow2(am2) + fem::pow2(pq21) + fem::pow2(pq22);
  pzcm = fem::sqrt(fem::abs(fem::pow2(amt) + fem::pow2(amt1) +
    fem::pow2(amt2) - 2.0f * amt * amt1 - 2.0f * amt * amt2 - 2.0f *
    amt1 * amt2)) / 2.0f / fem::sqrt(amt);
  //C                *******PZ of end-partons in c.m. frame of the string
  k(1, 1) = 2;
  k(1, 2) = kf1;
  p(1, 1) = pq11;
  p(1, 2) = pq12;
  p(1, 3) = pzcm;
  p(1, 4) = fem::sqrt(amt1 + fem::pow2(pzcm));
  p(1, 5) = am1;
  k(2, 1) = 1;
  k(2, 2) = kf2;
  p(2, 1) = pq21;
  p(2, 2) = pq22;
  p(2, 3) = -pzcm;
  p(2, 4) = fem::sqrt(amt2 + fem::pow2(pzcm));
  p(2, 5) = am2;
  n = 2;
  //C*****
  hirobo(0.0f, 0.0f, 0.0f, 0.0f, btz);
  jetot = 0;
  if ((fem::pow2(pq21) + fem::pow2(pq22)) > (fem::pow2(pq11) +
      fem::pow2(pq12))) {
    pmax1 = p(2, 1);
    pmax2 = p(2, 2);
    pmax3 = p(2, 3);
  }
  else {
    pmax1 = p(1, 1);
    pmax2 = p(1, 2);
    pmax3 = p(1, 3);
  }
  if (ntp == 1) {
    pp(jtp, 10) = pmax1;
    pp(jtp, 11) = pmax2;
    pp(jtp, 12) = pmax3;
  }
  else if (ntp == 2) {
    pt(jtp, 10) = pmax1;
    pt(jtp, 11) = pmax2;
    pt(jtp, 12) = pmax3;
  }
  //C*******************attach produced jets to the leadng partons****
  if (ntp == 1 && npj(jtp) != 0) {
    jetot = npj(jtp);
    //C                IF(NPJ(JTP).GE.2) CALL HIJSRT(JTP,1)
    //C                        ********sort jets in order of y
    iex = 0;
    if ((fem::abs(kf1) > 1000 && kf1 < 0) || (fem::abs(kf1) < 1000 &&
        kf1 > 0)) {
      iex = 1;
    }
    FEM_DOSTEP(i, n, 2, -1) {
      FEM_DO_SAFE(j, 1, 5) {
        ii = npj(jtp) + i;
        k(ii, j) = k(i, j);
        p(ii, j) = p(i, j);
        v(ii, j) = v(i, j);
      }
    }
    //C
    FEM_DO_SAFE(i, 1, npj(jtp)) {
      FEM_DO_SAFE(j, 1, 5) {
        k(i + 1, j) = 0;
        v(i + 1, j) = 0;
      }
      i0 = i;
      //Clin-4/12/01:                        IF(IEX.EQ.1) I0=NPJ(JTP)-I+1
      if (iex == 1 && (isoft != 2 || isflag != 0)) {
        i0 = npj(jtp) - i + 1;
      }
      //C                                ********reverse the order of jets
      kk1 = kfpj(jtp, i0);
      k(i + 1, 1) = 2;
      k(i + 1, 2) = kk1;
      if (kk1 != 21 && kk1 != 0) {
        k(i + 1, 1) = 1 + (fem::abs(kk1) + (2 * iex - 1) * kk1) / 2 /
          fem::abs(kk1);
      }
      p(i + 1, 1) = pjpx(jtp, i0);
      p(i + 1, 2) = pjpy(jtp, i0);
      p(i + 1, 3) = pjpz(jtp, i0);
      p(i + 1, 4) = pjpe(jtp, i0);
      p(i + 1, 5) = pjpm(jtp, i0);
    }
    n += npj(jtp);
  }
  else if (ntp == 2 && ntj(jtp) != 0) {
    jetot = ntj(jtp);
    //C                IF(NTJ(JTP).GE.2)  CALL HIJSRT(JTP,2)
    //C                        ********sort jets in order of y
    iex = 1;
    if ((fem::abs(kf2) > 1000 && kf2 < 0) || (fem::abs(kf2) < 1000 &&
        kf2 > 0)) {
      iex = 0;
    }
    FEM_DOSTEP(i, n, 2, -1) {
      FEM_DO_SAFE(j, 1, 5) {
        ii = ntj(jtp) + i;
        k(ii, j) = k(i, j);
        p(ii, j) = p(i, j);
        v(ii, j) = v(i, j);
      }
    }
    FEM_DO_SAFE(i, 1, ntj(jtp)) {
      FEM_DO_SAFE(j, 1, 5) {
        k(i + 1, j) = 0;
        v(i + 1, j) = 0;
      }
      i0 = i;
      //Clin-4/12/01:                        IF(IEX.EQ.1) I0=NTJ(JTP)-I+1
      if (iex == 1 && (isoft != 2 || isflag != 0)) {
        i0 = ntj(jtp) - i + 1;
      }
      //C                                ********reverse the order of jets
      kk1 = kftj(jtp, i0);
      k(i + 1, 1) = 2;
      k(i + 1, 2) = kk1;
      if (kk1 != 21 && kk1 != 0) {
        k(i + 1, 1) = 1 + (fem::abs(kk1) + (2 * iex - 1) * kk1) / 2 /
          fem::abs(kk1);
      }
      p(i + 1, 1) = pjtx(jtp, i0);
      p(i + 1, 2) = pjty(jtp, i0);
      p(i + 1, 3) = pjtz(jtp, i0);
      p(i + 1, 4) = pjte(jtp, i0);
      p(i + 1, 5) = pjtm(jtp, i0);
    }
    n += ntj(jtp);
  }
  if (ihpr2(1) > 0 && ranart(nseed) <= hidat(3)) {
    hdat20 = hidat(2);
    hpr150 = hipr1(5);
    if (ihpr2(8) == 0 && ihpr2(3) == 0 && ihpr2(9) == 0) {
      hidat(2) = 2.0f;
    }
    if (hint1(1) >= 1000.0f && jetot == 0) {
      hidat(2) = 3.0f;
      hipr1(5) = 5.0f;
    }
    attrad(cmn, ierror);
    hidat(2) = hdat20;
    hipr1(5) = hpr150;
  }
  else if (jetot == 0 && ihpr2(1) > 0 && hint1(1) >= 1000.0f &&
    ranart(nseed) <= 0.8f) {
    hdat20 = hidat(2);
    hpr150 = hipr1(5);
    hidat(2) = 3.0f;
    hipr1(5) = 5.0f;
    if (ihpr2(8) == 0 && ihpr2(3) == 0 && ihpr2(9) == 0) {
      hidat(2) = 2.0f;
    }
    attrad(cmn, ierror);
    hidat(2) = hdat20;
    hipr1(5) = hpr150;
  }
  if (ierror != 0) {
    return;
  }
  //C                ******** conduct soft radiations
  //C****************************
  //C
  //Clin-4/11/01 soft:
  //C        CALL LUEXEC
  if (isoft != 2 || isflag != 0) {
    luexec();
  }
  //C
  return;
  //C
  statement_1000:
  n = 1;
  k(1, 1) = 1;
  k(1, 2) = nfp(jtp, 3);
  FEM_DO_SAFE(jj, 1, 5) {
    p(1, jj) = pp(jtp, jj);
  }
  //C                        ********proj remain as a nucleon or delta
  //Clin-4/11/01 soft:
  //C        CALL LUEXEC
  if (isoft != 2 || isflag != 0) {
    luexec();
  }
  //C
  //C        call lulist(1)
  return;
  //C
  statement_1200:
  n = 1;
  k(1, 1) = 1;
  k(1, 2) = nft(jtp, 3);
  FEM_DO_SAFE(jj, 1, 5) {
    p(1, jj) = pt(jtp, jj);
  }
  //C                        ********targ remain as a nucleon or delta
  //Clin-4/11/01 soft:
  //C        CALL LUEXEC
  if (isoft != 2 || isflag != 0) {
    luexec();
  }
  //C
  //C        call lulist(1)
}

struct hijhrd_save
{
  float ampx;
  float amtx;
  float ecut1;
  float ecut2;
  float epm;
  float epp;
  float etm;
  float etp;
  int i;
  int iinird;
  int iopjet;
  arr<int, 2> ip;
  int ip1;
  int ip2;
  arr<int> ipb;
  arr<int> ipq;
  int is7;
  int is8;
  int isub11;
  int isub12;
  int isub28;
  arr<int, 2> it;
  int it1;
  int it2;
  arr<int> itb;
  arr<int> itq;
  int j;
  int jj;
  int jpp;
  int jtt;
  int l0;
  int lp;
  int lpb;
  int lpq;
  int lt;
  int ltb;
  int ltq;
  int misp;
  int miss;
  int mist;
  int mxjt;
  int mxsg;
  int mxsj;
  float pep;
  float pet;
  float pinird;
  float pxp;
  float pxt;
  float pyp;
  float pyt;
  float pzp;
  float pzt;
  float qm;
  float qmass2;
  float sw;
  float sxx;
  float wm;
  float wp;

  hijhrd_save() :
    ampx(fem::float0),
    amtx(fem::float0),
    ecut1(fem::float0),
    ecut2(fem::float0),
    epm(fem::float0),
    epp(fem::float0),
    etm(fem::float0),
    etp(fem::float0),
    i(fem::int0),
    iinird(fem::int0),
    iopjet(fem::int0),
    ip(dimension(100, 2), fem::fill0),
    ip1(fem::int0),
    ip2(fem::int0),
    ipb(dimension(50), fem::fill0),
    ipq(dimension(50), fem::fill0),
    is7(fem::int0),
    is8(fem::int0),
    isub11(fem::int0),
    isub12(fem::int0),
    isub28(fem::int0),
    it(dimension(100, 2), fem::fill0),
    it1(fem::int0),
    it2(fem::int0),
    itb(dimension(50), fem::fill0),
    itq(dimension(50), fem::fill0),
    j(fem::int0),
    jj(fem::int0),
    jpp(fem::int0),
    jtt(fem::int0),
    l0(fem::int0),
    lp(fem::int0),
    lpb(fem::int0),
    lpq(fem::int0),
    lt(fem::int0),
    ltb(fem::int0),
    ltq(fem::int0),
    misp(fem::int0),
    miss(fem::int0),
    mist(fem::int0),
    mxjt(fem::int0),
    mxsg(fem::int0),
    mxsj(fem::int0),
    pep(fem::float0),
    pet(fem::float0),
    pinird(fem::float0),
    pxp(fem::float0),
    pxt(fem::float0),
    pyp(fem::float0),
    pyt(fem::float0),
    pzp(fem::float0),
    pzt(fem::float0),
    qm(fem::float0),
    qmass2(fem::float0),
    sw(fem::float0),
    sxx(fem::float0),
    wm(fem::float0),
    wp(fem::float0)
  {}
};

void
hijhrd(
  common& cmn,
  int const& jp,
  int const& jt,
  int const& jout,
  int& jflg,
  int const& iopjt0)
{
  FEM_CMN_SVE(hijhrd);
  common_write write(cmn);
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  arr_ref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<float> hidat(cmn.hidat, dimension(10));
  arr_ref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_ref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  arr_ref<int> npj(cmn.npj, dimension(300));
  arr_ref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_ref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_ref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_ref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_ref<float, 2> pjpe(cmn.pjpe, dimension(300, 500));
  arr_ref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_ref<int> ntj(cmn.ntj, dimension(300));
  arr_ref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_ref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_ref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_ref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_ref<float, 2> pjte(cmn.pjte, dimension(300, 500));
  arr_ref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  int& nsg = cmn.nsg;
  const int maxstr = 150001;
  arr_ref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_ref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_ref<int, 2> k1sg(cmn.k1sg, dimension(maxstr, 100));
  arr_ref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_ref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_ref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_ref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_ref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_ref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  int& ndr = cmn.ndr;
  arr_ref<int, 2> iadr(cmn.iadr, dimension(maxstr, 2));
  arr_ref<int> kfdr(cmn.kfdr, dimension(maxstr));
  arr_ref<float, 2> pdr(cmn.pdr, dimension(maxstr, 5));
  arr_ref<float, 2> rtdr(cmn.rtdr, dimension(maxstr, 2));
  int& n = static_cast<common_lujets&>(cmn).n;
  arr_cref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_cref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_cref<float> vint(cmn.vint, dimension(400));
  arr_ref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  int& mint4 = cmn.mint4;
  int& mint5 = cmn.mint5;
  arr_cref<float, 2> atco(cmn.atco, dimension(200, 20));
  arr_cref<float> atxs(cmn.atxs, dim1(0, 200));
  float& pttrig = cmn.pttrig;
  int& maxmiss = cmn.maxmiss;
  //
  float& ampx = sve.ampx;
  float& amtx = sve.amtx;
  float& ecut1 = sve.ecut1;
  float& ecut2 = sve.ecut2;
  float& epm = sve.epm;
  float& epp = sve.epp;
  float& etm = sve.etm;
  float& etp = sve.etp;
  int& i = sve.i;
  int& iinird = sve.iinird;
  int& iopjet = sve.iopjet;
  arr_ref<int, 2> ip(sve.ip, dimension(100, 2));
  int& ip1 = sve.ip1;
  int& ip2 = sve.ip2;
  arr_ref<int> ipb(sve.ipb, dimension(50));
  arr_ref<int> ipq(sve.ipq, dimension(50));
  int& is7 = sve.is7;
  int& is8 = sve.is8;
  int& isub11 = sve.isub11;
  int& isub12 = sve.isub12;
  int& isub28 = sve.isub28;
  arr_ref<int, 2> it(sve.it, dimension(100, 2));
  int& it1 = sve.it1;
  int& it2 = sve.it2;
  arr_ref<int> itb(sve.itb, dimension(50));
  arr_ref<int> itq(sve.itq, dimension(50));
  int& j = sve.j;
  int& jj = sve.jj;
  int& jpp = sve.jpp;
  int& jtt = sve.jtt;
  int& l0 = sve.l0;
  int& lp = sve.lp;
  int& lpb = sve.lpb;
  int& lpq = sve.lpq;
  int& lt = sve.lt;
  int& ltb = sve.ltb;
  int& ltq = sve.ltq;
  int& misp = sve.misp;
  int& miss = sve.miss;
  int& mist = sve.mist;
  int& mxjt = sve.mxjt;
  int& mxsg = sve.mxsg;
  int& mxsj = sve.mxsj;
  float& pep = sve.pep;
  float& pet = sve.pet;
  float& pinird = sve.pinird;
  float& pxp = sve.pxp;
  float& pxt = sve.pxt;
  float& pyp = sve.pyp;
  float& pyt = sve.pyt;
  float& pzp = sve.pzp;
  float& pzt = sve.pzt;
  float& qm = sve.qm;
  float& qmass2 = sve.qmass2;
  float& sw = sve.sw;
  float& sxx = sve.sxx;
  float& wm = sve.wm;
  float& wp = sve.wp;
  //C
  //C        IOPTJET=1, ALL JET WILL FORM SINGLE STRING SYSTEM
  //C                0, ONLY Q-QBAR JET FORM SINGLE STRING SYSTEM
  //C*******Perform jets production and fragmentation when JP JT *******
  //C     scatter. JOUT-> number of hard scatterings precede this one  *
  //C     for the the same pair(JP,JT). JFLG->a flag to show whether   *
  //C     jets can be produced (with valence quark=1,gluon=2, q-qbar=3)*
  //C     or not(0). Information of jets are in  COMMON/ATTJET and     *
  //C     /MINJET. ABS(NFP(JP,6)) is the total number jets produced by *
  //C    JP. If NFP(JP,6)<0 JP can not produce jet anymore.                   *
  //C*******************************************************************
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HIJDAT/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /HJJET1/
  //Cc      SAVE /HJJET2/
  //C        COMMON/HJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
  //Cc      SAVE /HJJET4/
  //Cc      SAVE /RNDF77/
  //C************************************ HIJING common block
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /LUDAT1/
  //Cc      SAVE /PYSUBS/
  //Cc      SAVE /PYPARS/
  //Cc      SAVE /PYINT1/
  //Cc      SAVE /PYINT2/
  //Cc      SAVE /PYINT5/
  //Cc      SAVE /HPINT/
  //Clin-2/2012 correction:
  //C*********************************** LU common block
  mxjt = 500;
  //C                SIZE OF COMMON BLOCK FOR # OF PARTON PER STRING
  mxsg = 900;
  //C                SIZE OF COMMON BLOCK FOR # OF SINGLE STRINGS
  mxsj = 100;
  //C                SIZE OF COMMON BLOCK FOR # OF PARTON PER SINGLE
  //C                STRING
  jflg = 0;
  ihnt2(11) = jp;
  ihnt2(12) = jt;
  //C
  iopjet = iopjt0;
  if (iopjet == 1 && (nfp(jp, 6) != 0 || nft(jt, 6) != 0)) {
    iopjet = 0;
  }
  if (jp > ihnt2(1) || jt > ihnt2(3)) {
    return;
  }
  if (nfp(jp, 6) < 0 || nft(jt, 6) < 0) {
    return;
  }
  //C                ******** JP or JT can not produce jet anymore
  //C
  if (jout == 0) {
    epp = pp(jp, 4) + pp(jp, 3);
    epm = pp(jp, 4) - pp(jp, 3);
    etp = pt(jt, 4) + pt(jt, 3);
    etm = pt(jt, 4) - pt(jt, 3);
    if (epp < 0.0f) {
      goto statement_1000;
    }
    if (epm < 0.0f) {
      goto statement_1000;
    }
    if (etp < 0.0f) {
      goto statement_1000;
    }
    if (etm < 0.0f) {
      goto statement_1000;
    }
    if (epp / (epm + 0.01f) <= etp / (etm + 0.01f)) {
      return;
    }
  }
  //C                ********for the first hard scattering of (JP,JT)
  //C                        have collision only when Ycm(JP)>Ycm(JT)
  //C
  ecut1 = hipr1(1) + hipr1(8) + pp(jp, 14) + pp(jp, 15);
  ecut2 = hipr1(1) + hipr1(8) + pt(jt, 14) + pt(jt, 15);
  if (pp(jp, 4) <= ecut1) {
    nfp(jp, 6) = -fem::abs(nfp(jp, 6));
    return;
  }
  if (pt(jt, 4) <= ecut2) {
    nft(jt, 6) = -fem::abs(nft(jt, 6));
    return;
  }
  //C                *********must have enough energy to produce jets
  //C
  miss = 0;
  misp = 0;
  mist = 0;
  //C
  if (nfp(jp, 10) == 0 && nft(jt, 10) == 0) {
    mint(44) = mint4;
    mint(45) = mint5;
    xsec(0, 1) = atxs(0);
    xsec(11, 1) = atxs(11);
    xsec(12, 1) = atxs(12);
    xsec(28, 1) = atxs(28);
    FEM_DO_SAFE(i, 1, 20) {
      coef(11, i) = atco(11, i);
      coef(12, i) = atco(12, i);
      coef(28, i) = atco(28, i);
    }
  }
  else {
    isub11 = 0;
    isub12 = 0;
    isub28 = 0;
    if (xsec(11, 1) != 0) {
      isub11 = 1;
    }
    if (xsec(12, 1) != 0) {
      isub12 = 1;
    }
    if (xsec(28, 1) != 0) {
      isub28 = 1;
    }
    mint(44) = mint4 - isub11 - isub12 - isub28;
    mint(45) = mint5 - isub11 - isub12 - isub28;
    xsec(0, 1) = atxs(0) - atxs(11) - atxs(12) - atxs(28);
    xsec(11, 1) = 0.0f;
    xsec(12, 1) = 0.0f;
    xsec(28, 1) = 0.0f;
    FEM_DO_SAFE(i, 1, 20) {
      coef(11, i) = 0.0f;
      coef(12, i) = 0.0f;
      coef(28, i) = 0.0f;
    }
  }
  //C        ********Scatter the valence quarks only once per NN
  //C       collision,
  //C                afterwards only gluon can have hard scattering.
  statement_155:
  pythia();
  jj = mint(31);
  if (jj != 1) {
    goto statement_155;
  }
  //C                *********one hard collision at a time
  if (k(7, 2) ==  - k(8, 2)) {
    qmass2 = fem::pow2((p(7, 4) + p(8, 4))) - fem::pow2((p(7, 1) + p(8,
      1))) - fem::pow2((p(7, 2) + p(8, 2))) - fem::pow2((p(7, 3) + p(8,
      3)));
    qm = ulmass(k(7, 2));
    if (qmass2 < fem::pow2((2.0f * qm + hipr1(1)))) {
      goto statement_155;
    }
  }
  //C                ********q-qbar jets must has minimum mass HIPR1(1)
  pxp = pp(jp, 1) - p(3, 1);
  pyp = pp(jp, 2) - p(3, 2);
  pzp = pp(jp, 3) - p(3, 3);
  pep = pp(jp, 4) - p(3, 4);
  pxt = pt(jt, 1) - p(4, 1);
  pyt = pt(jt, 2) - p(4, 2);
  pzt = pt(jt, 3) - p(4, 3);
  pet = pt(jt, 4) - p(4, 4);
  //C
  if (pep <= ecut1) {
    misp++;
    if (misp < 50) {
      goto statement_155;
    }
    nfp(jp, 6) = -fem::abs(nfp(jp, 6));
    return;
  }
  if (pet <= ecut2) {
    mist++;
    if (mist < 50) {
      goto statement_155;
    }
    nft(jt, 6) = -fem::abs(nft(jt, 6));
    return;
  }
  //C                ******** if the remain energy<ECUT the proj or targ
  //C                         can not produce jet anymore
  //C
  wp = pep + pzp + pet + pzt;
  wm = pep - pzp + pet - pzt;
  if (wp < 0.0f || wm < 0.0f) {
    miss++;
    //Clin-6/2009 Let user set the limit when selecting high-Pt events
    //C     because more attempts may be needed:
    //C                IF(MISS.LT.50) GO TO 155
    if (pttrig > 0) {
      if (miss < maxmiss) {
        write(6, star), "Failed to generate minijet Pt>", pttrig, "GeV";
        goto statement_155;
      }
    }
    else {
      if (miss < 50) {
        goto statement_155;
      }
    }
    //C
    return;
  }
  //C                ********the total W+, W- must be positive
  sw = wp * wm;
  ampx = fem::sqrt(fem::pow2((ecut1 - hipr1(8))) + fem::pow2(pxp) +
    fem::pow2(pyp) + 0.01f);
  amtx = fem::sqrt(fem::pow2((ecut2 - hipr1(8))) + fem::pow2(pxt) +
    fem::pow2(pyt) + 0.01f);
  sxx = fem::pow2((ampx + amtx));
  if (sw < sxx || vint(43) < hipr1(1)) {
    miss++;
    //Clin-6/2009
    //C                IF(MISS.LT.50) GO TO 155
    if (miss > maxmiss) {
      goto statement_155;
    }
    return;
  }
  //C                ********the proj and targ remnants must have at least
  //C                        a CM energy that can produce two strings
  //C                        with minimum mass HIPR1(1)(see HIJSFT HIJFRG)
  //C
  hint1(41) = p(7, 1);
  hint1(42) = p(7, 2);
  hint1(43) = p(7, 3);
  hint1(44) = p(7, 4);
  hint1(45) = p(7, 5);
  hint1(46) = fem::sqrt(fem::pow2(p(7, 1)) + fem::pow2(p(7, 2)));
  hint1(51) = p(8, 1);
  hint1(52) = p(8, 2);
  hint1(53) = p(8, 3);
  hint1(54) = p(8, 4);
  hint1(55) = p(8, 5);
  hint1(56) = fem::sqrt(fem::pow2(p(8, 1)) + fem::pow2(p(8, 2)));
  ihnt2(14) = k(7, 2);
  ihnt2(15) = k(8, 2);
  //C
  pinird = (1.0f - fem::exp(-2.0f * (vint(47) - hidat(1)))) / (1.0f +
    fem::exp(-2.0f * (vint(47) - hidat(1))));
  iinird = 0;
  if (ranart(cmn.nseed) <= pinird) {
    iinird = 1;
  }
  if (k(7, 2) ==  - k(8, 2)) {
    goto statement_190;
  }
  if (k(7, 2) == 21 && k(8, 2) == 21 && iopjet == 1) {
    goto statement_190;
  }
  //C*******************************************************************
  //C        gluon  jets are going to be connectd with
  //C        the final leadng string of quark-aintquark
  //C*******************************************************************
  jflg = 2;
  jpp = 0;
  lpq = 0;
  lpb = 0;
  jtt = 0;
  ltq = 0;
  ltb = 0;
  is7 = 0;
  is8 = 0;
  hint1(47) = 0.0f;
  hint1(48) = 0.0f;
  hint1(49) = 0.0f;
  hint1(50) = 0.0f;
  hint1(67) = 0.0f;
  hint1(68) = 0.0f;
  hint1(69) = 0.0f;
  hint1(70) = 0.0f;
  FEM_DO_SAFE(i, 9, n) {
    if (k(i, 3) == 1 || k(i, 3) == 2 || fem::abs(k(i, 2)) > 30) {
      goto statement_180;
    }
    //C************************************************************
    if (k(i, 3) == 7) {
      hint1(47) += p(i, 1);
      hint1(48) += p(i, 2);
      hint1(49) += p(i, 3);
      hint1(50) += p(i, 4);
    }
    if (k(i, 3) == 8) {
      hint1(67) += p(i, 1);
      hint1(68) += p(i, 2);
      hint1(69) += p(i, 3);
      hint1(70) += p(i, 4);
    }
    //C************************modifcation made on Apr 10. 1996*****
    if (k(i, 2) > 21 && k(i, 2) <= 30) {
      ndr++;
      iadr(ndr, 1) = jp;
      iadr(ndr, 2) = jt;
      kfdr(ndr) = k(i, 2);
      pdr(ndr, 1) = p(i, 1);
      pdr(ndr, 2) = p(i, 2);
      pdr(ndr, 3) = p(i, 3);
      pdr(ndr, 4) = p(i, 4);
      pdr(ndr, 5) = p(i, 5);
      rtdr(ndr, 1) = 0.5f * (yp(1, jp) + yt(1, jt));
      rtdr(ndr, 2) = 0.5f * (yp(2, jp) + yt(2, jt));
      //C************************************************************
      goto statement_180;
      //C************************correction made on Oct. 14,1994*****
    }
    if (k(i, 3) == 7 || k(i, 3) == 3) {
      if (k(i, 3) == 7 && k(i, 2) != 21 && k(i, 2) == k(7, 2) && is7 == 0) {
        pp(jp, 10) = p(i, 1);
        pp(jp, 11) = p(i, 2);
        pp(jp, 12) = p(i, 3);
        pzp += p(i, 3);
        pep += p(i, 4);
        nfp(jp, 10) = 1;
        is7 = 1;
        goto statement_180;
      }
      if (k(i, 3) == 3 && (k(i, 2) != 21 || iinird == 0)) {
        pxp += p(i, 1);
        pyp += p(i, 2);
        pzp += p(i, 3);
        pep += p(i, 4);
        goto statement_180;
      }
      jpp++;
      ip(jpp, 1) = i;
      ip(jpp, 2) = 0;
      if (k(i, 2) != 21) {
        if (k(i, 2) > 0) {
          lpq++;
          ipq(lpq) = jpp;
          ip(jpp, 2) = lpq;
        }
        else if (k(i, 2) < 0) {
          lpb++;
          ipb(lpb) = jpp;
          ip(jpp, 2) = -lpb;
        }
      }
    }
    else if (k(i, 3) == 8 || k(i, 3) == 4) {
      if (k(i, 3) == 8 && k(i, 2) != 21 && k(i, 2) == k(8, 2) && is8 == 0) {
        pt(jt, 10) = p(i, 1);
        pt(jt, 11) = p(i, 2);
        pt(jt, 12) = p(i, 3);
        pzt += p(i, 3);
        pet += p(i, 4);
        nft(jt, 10) = 1;
        is8 = 1;
        goto statement_180;
      }
      if (k(i, 3) == 4 && (k(i, 2) != 21 || iinird == 0)) {
        pxt += p(i, 1);
        pyt += p(i, 2);
        pzt += p(i, 3);
        pet += p(i, 4);
        goto statement_180;
      }
      jtt++;
      it(jtt, 1) = i;
      it(jtt, 2) = 0;
      if (k(i, 2) != 21) {
        if (k(i, 2) > 0) {
          ltq++;
          itq(ltq) = jtt;
          it(jtt, 2) = ltq;
        }
        else if (k(i, 2) < 0) {
          ltb++;
          itb(ltb) = jtt;
          it(jtt, 2) = -ltb;
        }
      }
    }
    statement_180:;
  }
  //C
  if (lpq != lpb || ltq != ltb) {
    miss++;
    //Clin-6/2009
    //C                IF(MISS.LE.50) GO TO 155
    if (miss <= maxmiss) {
      goto statement_155;
    }
    write(6, star), " Q -QBAR NOT MATCHED IN HIJHRD";
    jflg = 0;
    return;
  }
  //C****The following will rearrange the partons so that a quark is***
  //C****allways followed by an anti-quark ****************************
  //C
  j = 0;
  statement_181:
  j++;
  if (j > jpp) {
    goto statement_182;
  }
  if (ip(j, 2) == 0) {
    goto statement_181;
  }
  else if (ip(j, 2) != 0) {
    lp = fem::abs(ip(j, 2));
    ip1 = ip(j, 1);
    ip2 = ip(j, 2);
    ip(j, 1) = ip(ipq(lp), 1);
    ip(j, 2) = ip(ipq(lp), 2);
    ip(ipq(lp), 1) = ip1;
    ip(ipq(lp), 2) = ip2;
    if (ip2 > 0) {
      ipq(ip2) = ipq(lp);
    }
    else if (ip2 < 0) {
      ipb(-ip2) = ipq(lp);
    }
    //C                ********replace J with a quark
    ip1 = ip(j + 1, 1);
    ip2 = ip(j + 1, 2);
    ip(j + 1, 1) = ip(ipb(lp), 1);
    ip(j + 1, 2) = ip(ipb(lp), 2);
    ip(ipb(lp), 1) = ip1;
    ip(ipb(lp), 2) = ip2;
    if (ip2 > 0) {
      ipq(ip2) = ipb(lp);
    }
    else if (ip2 < 0) {
      ipb(-ip2) = ipb(lp);
    }
    //C                ******** replace J+1 with anti-quark
    j++;
    goto statement_181;
  }
  //C
  statement_182:
  j = 0;
  statement_183:
  j++;
  if (j > jtt) {
    goto statement_184;
  }
  if (it(j, 2) == 0) {
    goto statement_183;
  }
  else if (it(j, 2) != 0) {
    lt = fem::abs(it(j, 2));
    it1 = it(j, 1);
    it2 = it(j, 2);
    it(j, 1) = it(itq(lt), 1);
    it(j, 2) = it(itq(lt), 2);
    it(itq(lt), 1) = it1;
    it(itq(lt), 2) = it2;
    if (it2 > 0) {
      itq(it2) = itq(lt);
    }
    else if (it2 < 0) {
      itb(-it2) = itq(lt);
    }
    //C                ********replace J with a quark
    it1 = it(j + 1, 1);
    it2 = it(j + 1, 2);
    it(j + 1, 1) = it(itb(lt), 1);
    it(j + 1, 2) = it(itb(lt), 2);
    it(itb(lt), 1) = it1;
    it(itb(lt), 2) = it2;
    if (it2 > 0) {
      itq(it2) = itb(lt);
    }
    else if (it2 < 0) {
      itb(-it2) = itb(lt);
    }
    //C                ******** replace J+1 with anti-quark
    j++;
    goto statement_183;
    //C
  }
  //C
  statement_184:
  if (npj(jp) + jpp > mxjt || ntj(jt) + jtt > mxjt) {
    jflg = 0;
    write(6, star), "number of partons per string exceeds";
    write(6, star), "the common block size";
    return;
  }
  //C                        ********check the bounds of common blocks
  FEM_DO_SAFE(j, 1, jpp) {
    kfpj(jp, npj(jp) + j) = k(ip(j, 1), 2);
    pjpx(jp, npj(jp) + j) = p(ip(j, 1), 1);
    pjpy(jp, npj(jp) + j) = p(ip(j, 1), 2);
    pjpz(jp, npj(jp) + j) = p(ip(j, 1), 3);
    pjpe(jp, npj(jp) + j) = p(ip(j, 1), 4);
    pjpm(jp, npj(jp) + j) = p(ip(j, 1), 5);
  }
  npj(jp) += jpp;
  FEM_DO_SAFE(j, 1, jtt) {
    kftj(jt, ntj(jt) + j) = k(it(j, 1), 2);
    pjtx(jt, ntj(jt) + j) = p(it(j, 1), 1);
    pjty(jt, ntj(jt) + j) = p(it(j, 1), 2);
    pjtz(jt, ntj(jt) + j) = p(it(j, 1), 3);
    pjte(jt, ntj(jt) + j) = p(it(j, 1), 4);
    pjtm(jt, ntj(jt) + j) = p(it(j, 1), 5);
  }
  ntj(jt) += jtt;
  goto statement_900;
  //C*****************************************************************
  //CThis is the case of a quark-antiquark jet it will fragment alone
  //C****************************************************************
  statement_190:
  jflg = 3;
  if (k(7, 2) != 21 && k(8, 2) != 21 && k(7, 2) * k(8, 2) > 0) {
    goto statement_155;
  }
  jpp = 0;
  lpq = 0;
  lpb = 0;
  FEM_DO_SAFE(i, 9, n) {
    if (k(i, 3) == 1 || k(i, 3) == 2 || fem::abs(k(i, 2)) > 30) {
      goto statement_200;
    }
    if (k(i, 2) > 21 && k(i, 2) <= 30) {
      ndr++;
      iadr(ndr, 1) = jp;
      iadr(ndr, 2) = jt;
      kfdr(ndr) = k(i, 2);
      pdr(ndr, 1) = p(i, 1);
      pdr(ndr, 2) = p(i, 2);
      pdr(ndr, 3) = p(i, 3);
      pdr(ndr, 4) = p(i, 4);
      pdr(ndr, 5) = p(i, 5);
      rtdr(ndr, 1) = 0.5f * (yp(1, jp) + yt(1, jt));
      rtdr(ndr, 2) = 0.5f * (yp(2, jp) + yt(2, jt));
      //C************************************************************
      goto statement_200;
      //C************************correction made on Oct. 14,1994*****
    }
    if (k(i, 3) == 3 && (k(i, 2) != 21 || iinird == 0)) {
      pxp += p(i, 1);
      pyp += p(i, 2);
      pzp += p(i, 3);
      pep += p(i, 4);
      goto statement_200;
    }
    if (k(i, 3) == 4 && (k(i, 2) != 21 || iinird == 0)) {
      pxt += p(i, 1);
      pyt += p(i, 2);
      pzt += p(i, 3);
      pet += p(i, 4);
      goto statement_200;
    }
    jpp++;
    ip(jpp, 1) = i;
    ip(jpp, 2) = 0;
    if (k(i, 2) != 21) {
      if (k(i, 2) > 0) {
        lpq++;
        ipq(lpq) = jpp;
        ip(jpp, 2) = lpq;
      }
      else if (k(i, 2) < 0) {
        lpb++;
        ipb(lpb) = jpp;
        ip(jpp, 2) = -lpb;
      }
    }
    statement_200:;
  }
  if (lpq != lpb) {
    miss++;
    //Clin-6/2009
    //C           IF(MISS.LE.50) GO TO 155
    if (miss <= maxmiss) {
      goto statement_155;
    }
    write(6, star), lpq, lpb, "Q-QBAR NOT CONSERVED OR NOT MATCHED";
    jflg = 0;
    return;
  }
  //C
  //C**** The following will rearrange the partons so that a quark is***
  //C**** allways followed by an anti-quark ****************************
  j = 0;
  statement_220:
  j++;
  if (j > jpp) {
    goto statement_222;
  }
  if (ip(j, 2) == 0) {
    goto statement_220;
  }
  lp = fem::abs(ip(j, 2));
  ip1 = ip(j, 1);
  ip2 = ip(j, 2);
  ip(j, 1) = ip(ipq(lp), 1);
  ip(j, 2) = ip(ipq(lp), 2);
  ip(ipq(lp), 1) = ip1;
  ip(ipq(lp), 2) = ip2;
  if (ip2 > 0) {
    ipq(ip2) = ipq(lp);
  }
  else if (ip2 < 0) {
    ipb(-ip2) = ipq(lp);
  }
  ipq(lp) = j;
  //C                ********replace J with a quark
  ip1 = ip(j + 1, 1);
  ip2 = ip(j + 1, 2);
  ip(j + 1, 1) = ip(ipb(lp), 1);
  ip(j + 1, 2) = ip(ipb(lp), 2);
  ip(ipb(lp), 1) = ip1;
  ip(ipb(lp), 2) = ip2;
  if (ip2 > 0) {
    ipq(ip2) = ipb(lp);
  }
  else if (ip2 < 0) {
    ipb(-ip2) = ipb(lp);
  }
  //C                ******** replace J+1 with an anti-quark
  ipb(lp) = j + 1;
  j++;
  goto statement_220;
  //C
  statement_222:
  if (lpq >= 1) {
    FEM_DO_SAFE(l0, 2, lpq) {
      ip1 = ip(2 * l0 - 3, 1);
      ip2 = ip(2 * l0 - 3, 2);
      ip(2 * l0 - 3, 1) = ip(ipq(l0), 1);
      ip(2 * l0 - 3, 2) = ip(ipq(l0), 2);
      ip(ipq(l0), 1) = ip1;
      ip(ipq(l0), 2) = ip2;
      if (ip2 > 0) {
        ipq(ip2) = ipq(l0);
      }
      else if (ip2 < 0) {
        ipb(-ip2) = ipq(l0);
      }
      ipq(l0) = 2 * l0 - 3;
      //C
      ip1 = ip(2 * l0 - 2, 1);
      ip2 = ip(2 * l0 - 2, 2);
      ip(2 * l0 - 2, 1) = ip(ipb(l0), 1);
      ip(2 * l0 - 2, 2) = ip(ipb(l0), 2);
      ip(ipb(l0), 1) = ip1;
      ip(ipb(l0), 2) = ip2;
      if (ip2 > 0) {
        ipq(ip2) = ipb(l0);
      }
      else if (ip2 < 0) {
        ipb(-ip2) = ipb(l0);
      }
      ipb(l0) = 2 * l0 - 2;
    }
    //C                ********move all the qqbar pair to the front of
    //C                                the list, except the first pair
    ip1 = ip(2 * lpq - 1, 1);
    ip2 = ip(2 * lpq - 1, 2);
    ip(2 * lpq - 1, 1) = ip(ipq(1), 1);
    ip(2 * lpq - 1, 2) = ip(ipq(1), 2);
    ip(ipq(1), 1) = ip1;
    ip(ipq(1), 2) = ip2;
    if (ip2 > 0) {
      ipq(ip2) = ipq(1);
    }
    else if (ip2 < 0) {
      ipb(-ip2) = ipq(1);
    }
    ipq(1) = 2 * lpq - 1;
    //C                ********move the first quark to the beginning of
    //C                                the last string system
    ip1 = ip(jpp, 1);
    ip2 = ip(jpp, 2);
    ip(jpp, 1) = ip(ipb(1), 1);
    ip(jpp, 2) = ip(ipb(1), 2);
    ip(ipb(1), 1) = ip1;
    ip(ipb(1), 2) = ip2;
    if (ip2 > 0) {
      ipq(ip2) = ipb(1);
    }
    else if (ip2 < 0) {
      ipb(-ip2) = ipb(1);
    }
    ipb(1) = jpp;
    //C                ********move the first anti-quark to the end of the
    //C                        last string system
  }
  if (nsg >= mxsg) {
    jflg = 0;
    write(6, star), "number of jets forming single strings exceeds";
    write(6, star), "the common block size";
    return;
  }
  if (jpp > mxsj) {
    jflg = 0;
    write(6, star), "number of partons per single jet system";
    write(6, star), "exceeds the common block size";
    return;
  }
  //C                ********check the bounds of common block size
  nsg++;
  njsg(nsg) = jpp;
  iasg(nsg, 1) = jp;
  iasg(nsg, 2) = jt;
  iasg(nsg, 3) = 0;
  FEM_DO_SAFE(i, 1, jpp) {
    k1sg(nsg, i) = 2;
    k2sg(nsg, i) = k(ip(i, 1), 2);
    if (k2sg(nsg, i) < 0) {
      k1sg(nsg, i) = 1;
    }
    pxsg(nsg, i) = p(ip(i, 1), 1);
    pysg(nsg, i) = p(ip(i, 1), 2);
    pzsg(nsg, i) = p(ip(i, 1), 3);
    pesg(nsg, i) = p(ip(i, 1), 4);
    pmsg(nsg, i) = p(ip(i, 1), 5);
  }
  k1sg(nsg, 1) = 2;
  k1sg(nsg, jpp) = 1;
  //C******* reset the energy-momentum of incoming particles ********
  statement_900:
  pp(jp, 1) = pxp;
  pp(jp, 2) = pyp;
  pp(jp, 3) = pzp;
  pp(jp, 4) = pep;
  pp(jp, 5) = 0.0f;
  pt(jt, 1) = pxt;
  pt(jt, 2) = pyt;
  pt(jt, 3) = pzt;
  pt(jt, 4) = pet;
  pt(jt, 5) = 0.0f;
  //C
  nfp(jp, 6)++;
  nft(jt, 6)++;
  return;
  //C
  statement_1000:
  jflg = -1;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), "Fatal HIJHRD error";
  write(6, star), jp, " proj E+,E-", epp, epm, " status", nfp(jp, 5);
  write(6, star), jt, " targ E+,E_", etp, etm, " status", nft(jt, 5);
}

struct jetini_save
{
  fem::str<16> beam;
  arr<float, 3> coef0;
  int i;
  int ilast;
  arr<int> ini;
  int isel;
  int isub;
  int itype;
  int j;
  arr<int> mint44;
  arr<int> mint45;
  fem::str<16> targ;
  arr<float, 2> xsec0;

  jetini_save() :
    beam(fem::char0),
    coef0(dimension(8, 200, 20), fem::fill0),
    i(fem::int0),
    ilast(fem::int0),
    ini(dimension(8), fem::fill0),
    isel(fem::int0),
    isub(fem::int0),
    itype(fem::int0),
    j(fem::int0),
    mint44(dimension(8), fem::fill0),
    mint45(dimension(8), fem::fill0),
    targ(fem::char0),
    xsec0(dim1(8).dim2(0, 200), fem::fill0)
  {}
};

void
jetini(
  common& cmn,
  int const& jp,
  int const& jt,
  int const& itrig)
{
  FEM_CMN_SVE(jetini);
  common_write write(cmn);
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_ref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_cref<int, 2> nft(cmn.nft, dimension(300, 15));
  int& mint4 = cmn.mint4;
  int& mint5 = cmn.mint5;
  arr_ref<float, 2> atco(cmn.atco, dimension(200, 20));
  arr_ref<float> atxs(cmn.atxs, dim1(0, 200));
  arr_cref<int, 2> mdcy(cmn.mdcy, dimension(500, 3));
  arr_ref<int, 2> mdme(cmn.mdme, dimension(2000, 2));
  int& msel = cmn.msel;
  arr_ref<int> msub(cmn.msub, dimension(200));
  arr_ref<float> ckin(cmn.ckin, dimension(200));
  arr_ref<int> mstp(cmn.mstp, dimension(200));
  arr_ref<float> parp(cmn.parp, dimension(200));
  arr_ref<int> mint(cmn.mint, dimension(400));
  arr_ref<float, 2> coef(cmn.coef, dimension(200, 20));
  arr_ref<float, 2> xsec(cmn.xsec, dim1(0, 200).dim2(3));
  //
  fem::str<16>& beam = sve.beam;
  arr_ref<float, 3> coef0(sve.coef0, dimension(8, 200, 20));
  int& i = sve.i;
  int& ilast = sve.ilast;
  arr_ref<int> ini(sve.ini, dimension(8));
  int& isel = sve.isel;
  int& isub = sve.isub;
  int& itype = sve.itype;
  int& j = sve.j;
  arr_ref<int> mint44(sve.mint44, dimension(8));
  arr_ref<int> mint45(sve.mint45, dimension(8));
  fem::str<16>& targ = sve.targ;
  arr_ref<float, 2> xsec0(sve.xsec0, dim1(8).dim2(0, 200));
  if (is_called_first_time) {
    fem::data((values, 8*datum(0))), ini;
    ilast = -1;
  }
  //C*******Initialize PYTHIA for jet production**********************
  //C        itrig=0: for normal processes
  //C        itrig=1: for triggered processes
  //C       JP: sequence number of the projectile
  //C       JT: sequence number of the target
  //C     For A+A collisions, one has to initilize pythia
  //C     separately for each type of collisions, pp, pn,np and nn,
  //C     or hp and hn for hA collisions. In this subroutine we use the following
  //C     catalogue for different type of collisions:
  //C     h+h: h+h (itype=1)
  //C     h+A: h+p (itype=1), h+n (itype=2)
  //C     A+h: p+h (itype=1), n+h (itype=2)
  //C     A+A: p+p (itype=1), p+n (itype=2), n+p (itype=3), n+n (itype=4)
  //C*****************************************************************
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /HPINT/
  //C
  //Cc      SAVE /LUDAT1/
  //Cc      SAVE /LUDAT3/
  //Cc      SAVE /PYSUBS/
  //Cc      SAVE /PYPARS/
  //Cc      SAVE /PYINT1/
  //Cc      SAVE /PYINT2/
  //Cc      SAVE /PYINT5/
  //Clin        DATA INI/8*0/ilast/-1/
  //C
  ihnt2(11) = jp;
  ihnt2(12) = jt;
  if (ihnt2(5) != 0 && ihnt2(6) != 0) {
    itype = 1;
  }
  else if (ihnt2(5) != 0 && ihnt2(6) == 0) {
    itype = 1;
    if (nft(jt, 4) == 2112) {
      itype = 2;
    }
  }
  else if (ihnt2(5) == 0 && ihnt2(6) != 0) {
    itype = 1;
    if (nfp(jp, 4) == 2112) {
      itype = 2;
    }
  }
  else {
    if (nfp(jp, 4) == 2212 && nft(jt, 4) == 2212) {
      itype = 1;
    }
    else if (nfp(jp, 4) == 2212 && nft(jt, 4) == 2112) {
      itype = 2;
    }
    else if (nfp(jp, 4) == 2112 && nft(jt, 4) == 2212) {
      itype = 3;
    }
    else {
      itype = 4;
    }
  }
  //C
  //Clin-12/2012 correct NN differential cross section in HIJING:
  //C        write(94,*) 'In JETINI: ',jp,jt,NFP(JP,4),NFT(JT,4),itype
  //C
  if (itrig != 0) {
    goto statement_160;
  }
  if (itrig == ilast) {
    goto statement_150;
  }
  mstp(2) = 2;
  //C                        ********second order running alpha_strong
  mstp(33) = 1;
  parp(31) = hipr1(17);
  //C                        ********inclusion of K factor
  mstp(51) = 3;
  //C                        ********Duke-Owens set 1 structure functions
  mstp(61) = 1;
  //C                        ********INITIAL STATE RADIATION
  mstp(71) = 1;
  //C                        ********FINAL STATE RADIATION
  if (ihpr2(2) == 0 || ihpr2(2) == 2) {
    mstp(61) = 0;
  }
  if (ihpr2(2) == 0 || ihpr2(2) == 1) {
    mstp(71) = 0;
  }
  //C
  mstp(81) = 0;
  //C                        ******** NO MULTIPLE INTERACTION
  mstp(82) = 1;
  //C                        *******STRUCTURE OF MUTLIPLE INTERACTION
  mstp(111) = 0;
  //C                ********frag off(have to be done by local call)
  if (ihpr2(10) == 0) {
    mstp(122) = 0;
  }
  //C                ********No printout of initialization information
  parp(81) = hipr1(8);
  ckin(5) = hipr1(8);
  ckin(3) = hipr1(8);
  ckin(4) = hipr1(9);
  if (hipr1(9) <= hipr1(8)) {
    ckin(4) = -1.0f;
  }
  ckin(9) = -10.0f;
  ckin(10) = 10.0f;
  msel = 0;
  FEM_DO_SAFE(isub, 1, 200) {
    msub(isub) = 0;
  }
  msub(11) = 1;
  msub(12) = 1;
  msub(13) = 1;
  msub(28) = 1;
  msub(53) = 1;
  msub(68) = 1;
  msub(81) = 1;
  msub(82) = 1;
  FEM_DO_SAFE(j, 1, fem::min(8, mdcy(21, 3))) {
    mdme(mdcy(21, 2) + j - 1, 1) = 0;
  }
  isel = 4;
  if (hint1(1) >= 20.0f && ihpr2(18) == 1) {
    isel = 5;
  }
  mdme(mdcy(21, 2) + isel - 1, 1) = 1;
  //C                        ********QCD subprocesses
  msub(14) = 1;
  msub(18) = 1;
  msub(29) = 1;
  //C                       ******* direct photon production
  statement_150:
  if (ini(itype) != 0) {
    goto statement_800;
  }
  goto statement_400;
  //C
  //C        *****triggered subprocesses, jet, photon, heavy quark and DY
  //C
  statement_160:
  itype += 4;
  if (itrig == ilast) {
    goto statement_260;
  }
  parp(81) = fem::abs(hipr1(10)) - 0.25f;
  ckin(5) = fem::abs(hipr1(10)) - 0.25f;
  ckin(3) = fem::abs(hipr1(10)) - 0.25f;
  ckin(4) = fem::abs(hipr1(10)) + 0.25f;
  if (hipr1(10) < hipr1(8)) {
    ckin(4) = -1.0f;
  }
  //C
  msel = 0;
  FEM_DO_SAFE(isub, 1, 200) {
    msub(isub) = 0;
  }
  if (ihpr2(3) == 1) {
    msub(11) = 1;
    msub(12) = 1;
    msub(13) = 1;
    msub(28) = 1;
    msub(53) = 1;
    msub(68) = 1;
    msub(81) = 1;
    msub(82) = 1;
    msub(14) = 1;
    msub(18) = 1;
    msub(29) = 1;
    FEM_DO_SAFE(j, 1, fem::min(8, mdcy(21, 3))) {
      mdme(mdcy(21, 2) + j - 1, 1) = 0;
    }
    isel = 4;
    if (hint1(1) >= 20.0f && ihpr2(18) == 1) {
      isel = 5;
    }
    mdme(mdcy(21, 2) + isel - 1, 1) = 1;
    //C                        ********QCD subprocesses
  }
  else if (ihpr2(3) == 2) {
    msub(14) = 1;
    msub(18) = 1;
    msub(29) = 1;
    //C                ********Direct photon production
    //C                q+qbar->g+gamma,q+qbar->gamma+gamma, q+g->q+gamma
  }
  else if (ihpr2(3) == 3) {
    ckin(3) = fem::max(0.0f, hipr1(10));
    ckin(5) = hipr1(8);
    parp(81) = hipr1(8);
    msub(81) = 1;
    msub(82) = 1;
    FEM_DO_SAFE(j, 1, fem::min(8, mdcy(21, 3))) {
      mdme(mdcy(21, 2) + j - 1, 1) = 0;
    }
    isel = 4;
    if (hint1(1) >= 20.0f && ihpr2(18) == 1) {
      isel = 5;
    }
    mdme(mdcy(21, 2) + isel - 1, 1) = 1;
    //C             **********Heavy quark production
  }
  statement_260:
  if (ini(itype) != 0) {
    goto statement_800;
  }
  //C
  statement_400:
  ini(itype) = 1;
  if (ihpr2(10) == 0) {
    mstp(122) = 0;
  }
  if (nfp(jp, 4) == 2212) {
    beam = "P";
  }
  else if (nfp(jp, 4) ==  - 2212) {
    beam = "P~";
  }
  else if (nfp(jp, 4) == 2112) {
    beam = "N";
  }
  else if (nfp(jp, 4) ==  - 2112) {
    beam = "N~";
  }
  else if (nfp(jp, 4) == 211) {
    beam = "PI+";
  }
  else if (nfp(jp, 4) ==  - 211) {
    beam = "PI-";
  }
  else if (nfp(jp, 4) == 321) {
    beam = "PI+";
  }
  else if (nfp(jp, 4) ==  - 321) {
    beam = "PI-";
  }
  else {
    write(6, star), "unavailable beam type", nfp(jp, 4);
  }
  if (nft(jt, 4) == 2212) {
    targ = "P";
  }
  else if (nft(jt, 4) ==  - 2212) {
    targ = "P~";
  }
  else if (nft(jt, 4) == 2112) {
    targ = "N";
  }
  else if (nft(jt, 4) ==  - 2112) {
    targ = "N~";
  }
  else if (nft(jt, 4) == 211) {
    targ = "PI+";
  }
  else if (nft(jt, 4) ==  - 211) {
    targ = "PI-";
  }
  else if (nft(jt, 4) == 321) {
    targ = "PI+";
  }
  else if (nft(jt, 4) ==  - 321) {
    targ = "PI-";
  }
  else {
    write(6, star), "unavailable target type", nft(jt, 4);
  }
  //C
  ihnt2(16) = 1;
  //C       ******************indicate for initialization use when
  //C                         structure functions are called in PYTHIA
  //C
  pyinit("CMS", beam, targ, hint1(1));
  mint4 = mint(44);
  mint5 = mint(45);
  mint44(itype) = mint(44);
  mint45(itype) = mint(45);
  atxs(0) = xsec(0, 1);
  xsec0(itype, 0) = xsec(0, 1);
  FEM_DO_SAFE(i, 1, 200) {
    atxs(i) = xsec(i, 1);
    xsec0(itype, i) = xsec(i, 1);
    FEM_DO_SAFE(j, 1, 20) {
      atco(i, j) = coef(i, j);
      coef0(itype, i, j) = coef(i, j);
    }
  }
  //C
  ihnt2(16) = 0;
  //C
  return;
  //C                ********Store the initialization information for
  //C                                late use
  //C
  statement_800:
  mint(44) = mint44(itype);
  mint(45) = mint45(itype);
  mint4 = mint(44);
  mint5 = mint(45);
  xsec(0, 1) = xsec0(itype, 0);
  atxs(0) = xsec(0, 1);
  FEM_DO_SAFE(i, 1, 200) {
    xsec(i, 1) = xsec0(itype, i);
    atxs(i) = xsec(i, 1);
    FEM_DO_SAFE(j, 1, 20) {
      coef(i, j) = coef0(itype, i, j);
      atco(i, j) = coef(i, j);
    }
  }
  ilast = itrig;
  mint(11) = nfp(jp, 4);
  mint(12) = nft(jt, 4);
}

struct attflv_save
{
  int id0;
  int id00;
  int nsign;
  float x;

  attflv_save() :
    id0(fem::int0),
    id00(fem::int0),
    nsign(fem::int0),
    x(fem::float0)
  {}
};

void
attflv(
  common& cmn,
  int const& id,
  int& idq,
  int& idqq)
{
  FEM_CMN_SVE(attflv);
  int& id0 = sve.id0;
  int& id00 = sve.id00;
  int& nsign = sve.nsign;
  float& x = sve.x;
  //Cc      SAVE /RNDF77/
  //C
  if (fem::abs(id) < 100) {
    nsign = 1;
    idq = id / 100;
    idqq = -id / 10 + idq * 10;
    if (fem::abs(idq) == 3) {
      nsign = -1;
    }
    idq = nsign * idq;
    idqq = nsign * idqq;
    if (idq < 0) {
      id0 = idq;
      idq = idqq;
      idqq = id0;
    }
    return;
  }
  //C                ********return ID of quark(IDQ) and anti-quark(IDQQ)
  //C                        for pions and kaons
  //C
  //C        Return LU ID for quarks and diquarks for proton(ID=2212)
  //C        anti-proton(ID=-2212) and nuetron(ID=2112)
  //C        LU ID for d=1,u=2, (ud)0=2101, (ud)1=2103,
  //C       (dd)1=1103,(uu)1=2203.
  //C        Use SU(6)  weight  proton=1/3d(uu)1 + 1/6u(ud)1 + 1/2u(ud)0
  //C                          nurtron=1/3u(dd)1 + 1/6d(ud)1 + 1/2d(ud)0
  //C
  idq = 2;
  if (fem::abs(id) == 2112) {
    idq = 1;
  }
  idqq = 2101;
  x = ranart(cmn.nseed);
  if (x <= 0.5f) {
    goto statement_30;
  }
  if (x > 0.666667f) {
    goto statement_10;
  }
  idqq = 2103;
  goto statement_30;
  statement_10:
  idq = 1;
  idqq = 2203;
  if (fem::abs(id) == 2112) {
    idq = 2;
    idqq = 1103;
  }
  statement_30:
  if (id < 0) {
    id00 = idqq;
    idqq = -idq;
    idq = -id00;
  }
}

struct hijini_save
{
  int i;
  int idq;
  int idqq;
  int ipp;
  int ipt;

  hijini_save() :
    i(fem::int0),
    idq(fem::int0),
    idqq(fem::int0),
    ipp(fem::int0),
    ipt(fem::int0)
  {}
};

void
hijini(
  common& cmn)
{
  FEM_CMN_SVE(hijini);
  // COMMON hparnt
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  // COMMON hstrng
  arr_ref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_ref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  // COMMON hjjet1
  arr_ref<int> npj(cmn.npj, dimension(300));
  arr_ref<int> ntj(cmn.ntj, dimension(300));
  // COMMON rndf77
  int& nseed = cmn.nseed;
  //
  // SAVE
  int& i = sve.i;
  int& idq = sve.idq;
  int& idqq = sve.idqq;
  int& ipp = sve.ipp;
  int& ipt = sve.ipt;
  //
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /HJJET1/
  //Cc      SAVE /HJJET2/
  //C        COMMON/HJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
  //Cc      SAVE /HJJET4/
  //Cc      SAVE /RNDF77/
  //C****************Reset the momentum of initial particles************
  //C             and assign flavors to the proj and targ string       *
  //C*******************************************************************
  cmn.nsg = 0;
  cmn.ndr = 0;
  ipp = 2212;
  ipt = 2212;
  if (ihnt2(5) != 0) {
    ipp = ihnt2(5);
  }
  if (ihnt2(6) != 0) {
    ipt = ihnt2(6);
  }
  //C                ********in case the proj or targ is a hadron.
  //C
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    pp(i, 1) = 0.0f;
    pp(i, 2) = 0.0f;
    pp(i, 3) = fem::sqrt(fem::pow2(hint1(1)) / 4.0f - fem::pow2(hint1(8)));
    pp(i, 4) = hint1(1) / 2;
    pp(i, 5) = hint1(8);
    pp(i, 6) = 0.0f;
    pp(i, 7) = 0.0f;
    pp(i, 8) = 0.0f;
    pp(i, 9) = 0.0f;
    pp(i, 10) = 0.0f;
    //Cbzdbg2/22/99
    //Ctest OFF
    pp(i, 11) = 0.0f;
    pp(i, 12) = 0.0f;
    //Cbzdbg2/22/99end
    nfp(i, 3) = ipp;
    nfp(i, 4) = ipp;
    nfp(i, 5) = 0;
    nfp(i, 6) = 0;
    nfp(i, 7) = 0;
    nfp(i, 8) = 0;
    nfp(i, 9) = 0;
    nfp(i, 10) = 0;
    nfp(i, 11) = 0;
    npj(i) = 0;
    if (i > fem::abs(ihnt2(2))) {
      nfp(i, 3) = 2112;
    }
    //C
    //Clin-12/2012 correct NN differential cross section in HIJING:
    if (i > fem::abs(ihnt2(2))) {
      nfp(i, 4) = 2112;
    }
    //C
    attflv(cmn, nfp(i, 3), idq, idqq);
    nfp(i, 1) = idq;
    nfp(i, 2) = idqq;
    nfp(i, 15) = -1;
    if (fem::abs(idq) > 1000 || (fem::abs(idq * idqq) < 100 && ranart(
        nseed) < 0.5f)) {
      nfp(i, 15) = 1;
    }
    pp(i, 14) = ulmass(idq);
    pp(i, 15) = ulmass(idqq);
  }
  //C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    pt(i, 1) = 0.0f;
    pt(i, 2) = 0.0f;
    pt(i, 3) = -fem::sqrt(fem::pow2(hint1(1)) / 4.0f - fem::pow2(hint1(9)));
    pt(i, 4) = hint1(1) / 2.0f;
    pt(i, 5) = hint1(9);
    pt(i, 6) = 0.0f;
    pt(i, 7) = 0.0f;
    pt(i, 8) = 0.0f;
    pt(i, 9) = 0.0f;
    pt(i, 10) = 0.0f;
    //Ctest OFF
    //Cbzdbg2/22/99
    pt(i, 11) = 0.0f;
    pt(i, 12) = 0.0f;
    //Cbzdbg2/22/99end
    nft(i, 3) = ipt;
    nft(i, 4) = ipt;
    nft(i, 5) = 0;
    nft(i, 6) = 0;
    nft(i, 7) = 0;
    nft(i, 8) = 0;
    nft(i, 9) = 0;
    nft(i, 10) = 0;
    nft(i, 11) = 0;
    ntj(i) = 0;
    if (i > fem::abs(ihnt2(4))) {
      nft(i, 3) = 2112;
    }
    //C
    //Clin-12/2012 correct NN differential cross section in HIJING:
    if (i > fem::abs(ihnt2(4))) {
      nft(i, 4) = 2112;
    }
    //C
    attflv(cmn, nft(i, 3), idq, idqq);
    nft(i, 1) = idq;
    nft(i, 2) = idqq;
    nft(i, 15) = 1;
    if (fem::abs(idq) > 1000 || (fem::abs(idq * idqq) < 100 && ranart(
        nseed) < 0.5f)) {
      nft(i, 15) = -1;
    }
    pt(i, 14) = ulmass(idq);
    pt(i, 15) = ulmass(idqq);
  }
}

struct hijels_save
{
  float am1;
  float am2;
  float amm;
  float bb;
  float cc;
  double db;
  double dbp;
  double dbx;
  double dby;
  double dbz;
  double dga;
  double dgabp;
  double dp1;
  double dp2;
  double dp3;
  double dp4;
  float ecm;
  float els;
  float els0;
  float ep;
  float pcm1;
  float pcm2;
  float pcm3;
  float phi;
  float pmax;
  float rr;
  float tt;

  hijels_save() :
    am1(fem::float0),
    am2(fem::float0),
    amm(fem::float0),
    bb(fem::float0),
    cc(fem::float0),
    db(fem::double0),
    dbp(fem::double0),
    dbx(fem::double0),
    dby(fem::double0),
    dbz(fem::double0),
    dga(fem::double0),
    dgabp(fem::double0),
    dp1(fem::double0),
    dp2(fem::double0),
    dp3(fem::double0),
    dp4(fem::double0),
    ecm(fem::float0),
    els(fem::float0),
    els0(fem::float0),
    ep(fem::float0),
    pcm1(fem::float0),
    pcm2(fem::float0),
    pcm3(fem::float0),
    phi(fem::float0),
    pmax(fem::float0),
    rr(fem::float0),
    tt(fem::float0)
  {}
};

//C
//C*******************************************************************
//CThis subroutine performs elastic scattering between two nucleons
//C
//C*******************************************************************
void
hijels(
  common& cmn,
  arr_ref<float> psc1,
  arr_ref<float> psc2)
{
  FEM_CMN_SVE(hijels);
  psc1(dimension(5));
  psc2(dimension(5));
  common_write write(cmn);
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  int& nseed = cmn.nseed;
  //
  float& am1 = sve.am1;
  float& am2 = sve.am2;
  float& amm = sve.amm;
  float& bb = sve.bb;
  float& cc = sve.cc;
  double& db = sve.db;
  double& dbp = sve.dbp;
  double& dbx = sve.dbx;
  double& dby = sve.dby;
  double& dbz = sve.dbz;
  double& dga = sve.dga;
  double& dgabp = sve.dgabp;
  double& dp1 = sve.dp1;
  double& dp2 = sve.dp2;
  double& dp3 = sve.dp3;
  double& dp4 = sve.dp4;
  float& ecm = sve.ecm;
  float& els = sve.els;
  float& els0 = sve.els0;
  float& ep = sve.ep;
  float& pcm1 = sve.pcm1;
  float& pcm2 = sve.pcm2;
  float& pcm3 = sve.pcm3;
  float& phi = sve.phi;
  float& pmax = sve.pmax;
  float& rr = sve.rr;
  float& tt = sve.tt;
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /RNDF77/
  //C
  cc = 1.0f - hint1(12) / hint1(13);
  rr = (1.0f - cc) * hint1(13) / hint1(12) / (1.0f - hipr1(33)) - 1.0f;
  bb = 0.5f * (3.0f + rr + fem::sqrt(9.0f + 10.0f * rr + fem::pow2(rr)));
  ep = fem::sqrt(fem::pow2((psc1(1) - psc2(1))) + fem::pow2((psc1(
    2) - psc2(2))) + fem::pow2((psc1(3) - psc2(3))));
  if (ep <= 0.1f) {
    return;
  }
  els0 = 98.0f / ep + 52.0f * fem::pow2((1.0f + rr));
  pcm1 = psc1(1) + psc2(1);
  pcm2 = psc1(2) + psc2(2);
  pcm3 = psc1(3) + psc2(3);
  ecm = psc1(4) + psc2(4);
  am1 = fem::pow2(psc1(5));
  am2 = fem::pow2(psc2(5));
  amm = fem::pow2(ecm) - fem::pow2(pcm1) - fem::pow2(pcm2) - fem::pow2(pcm3);
  if (amm <= psc1(5) + psc2(5)) {
    return;
  }
  //C                ********elastic scattering only when approaching
  //C                                to each other
  pmax = (fem::pow2(amm) + fem::pow2(am1) + fem::pow2(am2) - 2.0f *
    amm * am1 - 2.0f * amm * am2 - 2.0f * am1 * am2) / 4.0f / amm;
  pmax = fem::abs(pmax);
  statement_20:
  tt = ranart(nseed) * fem::min(pmax, 1.5f);
  els = 98.0f * fem::exp(-2.8f * tt) / ep + 52.0f * fem::exp(-9.2f *
    tt) * fem::pow2((1.0f + rr * fem::exp(-4.6f * (bb - 1.0f) * tt)));
  if (ranart(nseed) > els / els0) {
    goto statement_20;
  }
  phi = 2.0f * hipr1(40) * ranart(nseed);
  //C
  dbx = fem::dble(pcm1 / ecm);
  dby = fem::dble(pcm2 / ecm);
  dbz = fem::dble(pcm3 / ecm);
  db = fem::dsqrt(fem::pow2(dbx) + fem::pow2(dby) + fem::pow2(dbz));
  if (db > 0.99999999e0) {
    dbx = dbx * (0.99999999e0 / db);
    dby = dby * (0.99999999e0 / db);
    dbz = dbz * (0.99999999e0 / db);
    db = 0.99999999e0;
    write(6, star), " (HIJELS) boost vector too large";
    //C                ********Rescale boost vector if too close to unity.
  }
  dga = 1e0 / fem::sqrt(1e0 - fem::pow2(db));
  //C
  dp1 = fem::dble(fem::sqrt(tt) * fem::sin(phi));
  dp2 = fem::dble(fem::sqrt(tt) * fem::cos(phi));
  dp3 = fem::dble(fem::sqrt(pmax - tt));
  dp4 = fem::dble(fem::sqrt(pmax + am1));
  dbp = dbx * dp1 + dby * dp2 + dbz * dp3;
  dgabp = dga * (dga * dbp / (1e0 + dga) + dp4);
  psc1(1) = fem::sngl(dp1 + dgabp * dbx);
  psc1(2) = fem::sngl(dp2 + dgabp * dby);
  psc1(3) = fem::sngl(dp3 + dgabp * dbz);
  psc1(4) = fem::sngl(dga * (dp4 + dbp));
  //C
  dp1 = -fem::dble(fem::sqrt(tt) * fem::sin(phi));
  dp2 = -fem::dble(fem::sqrt(tt) * fem::cos(phi));
  dp3 = -fem::dble(fem::sqrt(pmax - tt));
  dp4 = fem::dble(fem::sqrt(pmax + am2));
  dbp = dbx * dp1 + dby * dp2 + dbz * dp3;
  dgabp = dga * (dga * dbp / (1e0 + dga) + dp4);
  psc2(1) = fem::sngl(dp1 + dgabp * dbx);
  psc2(2) = fem::sngl(dp2 + dgabp * dby);
  psc2(3) = fem::sngl(dp3 + dgabp * dbz);
  psc2(4) = fem::sngl(dga * (dp4 + dbp));
}

typedef float (*bk_function_pointer)(common&, float const&);

float
bk(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  //Cc      SAVE /BESEL/
  return_value = fem::exp(-x) * fem::pow((fem::pow2(x) - fem::pow2(cmn.x4)),
    2.50f) / 15.0f;
  return return_value;
}

struct gauss2_save
{
  float aa;
  float bb;
  float c1;
  float c2;
  float identifier_const;
  float delta;
  int i;
  float s16;
  float s8;
  float u;
  arr<float> w;
  arr<float> x;
  float y;

  gauss2_save() :
    aa(fem::float0),
    bb(fem::float0),
    c1(fem::float0),
    c2(fem::float0),
    identifier_const(fem::float0),
    delta(fem::float0),
    i(fem::int0),
    s16(fem::float0),
    s8(fem::float0),
    u(fem::float0),
    w(dimension(12), fem::fill0),
    x(dimension(12), fem::fill0),
    y(fem::float0)
  {}
};

float
gauss2(
  common& cmn,
  bk_function_pointer f,
  float const& a,
  float const& b,
  float const& eps)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(gauss2);
  common_write write(cmn);
  float& aa = sve.aa;
  float& bb = sve.bb;
  float& c1 = sve.c1;
  float& c2 = sve.c2;
  float& identifier_const = sve.identifier_const;
  float& delta = sve.delta;
  int& i = sve.i;
  float& s16 = sve.s16;
  float& s8 = sve.s8;
  float& u = sve.u;
  arr_ref<float> w(sve.w, dimension(12));
  arr_ref<float> x(sve.x, dimension(12));
  float& y = sve.y;
  if (is_called_first_time) {
    identifier_const = 1.0e-12f;
    {
      static const float values[] = {
        0.1012285f, .2223810f, .3137067f, .3623838f, .0271525f,
          .0622535f, 0.0951585f, .1246290f, .1495960f, .1691565f,
          .1826034f, .1894506f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        w;
    }
    {
      static const float values[] = {
        0.9602899f, .7966665f, .5255324f, .1834346f, .9894009f,
          .9445750f, 0.8656312f, .7554044f, .6178762f, .4580168f,
          .2816036f, .0950125f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        x;
    }
  }
  //C
  delta = identifier_const * fem::abs(a - b);
  return_value = 0.0f;
  aa = a;
  statement_5:
  y = b - aa;
  if (fem::abs(y) <= delta) {
    return return_value;
  }
  statement_2:
  bb = aa + y;
  c1 = 0.5f * (aa + bb);
  c2 = c1 - aa;
  s8 = 0.0f;
  s16 = 0.0f;
  FEM_DO_SAFE(i, 1, 4) {
    u = x(i) * c2;
    s8 += w(i) * (f(cmn, c1 + u) + f(cmn, c1 - u));
  }
  FEM_DO_SAFE(i, 5, 12) {
    u = x(i) * c2;
    s16 += w(i) * (f(cmn, c1 + u) + f(cmn, c1 - u));
  }
  s8 = s8 * c2;
  s16 = s16 * c2;
  if (fem::abs(s16 - s8) > eps * (1.f + fem::abs(s16))) {
    goto statement_4;
  }
  return_value += s16;
  aa = bb;
  goto statement_5;
  statement_4:
  y = 0.5f * y;
  if (fem::abs(y) > delta) {
    goto statement_2;
  }
  write(6, "(1x,'GAUSS2....TOO HIGH ACURACY REQUIRED')");
  return_value = 0.0f;
  return return_value;
}

float
omg0(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  // COMMON besel
  float& x4 = cmn.x4;
  //
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /BESEL/
  x4 = hipr1(32) * fem::sqrt(x);
  return_value = fem::pow2(hipr1(32)) * gauss2(cmn, bk, x4, x4 + 20.0f,
    0.01f) / 96.0f;
  return return_value;
}

struct romg_save
{
  arr<float> fr;
  int i;
  int i0;
  int ix;
  float xr;

  romg_save() :
    fr(dim1(0, 1000), fem::fill0),
    i(fem::int0),
    i0(fem::int0),
    ix(fem::int0),
    xr(fem::float0)
  {}
};

float
romg(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(romg);
  arr_ref<float> fr(sve.fr, dim1(0, 1000));
  int& i = sve.i;
  int& i0 = sve.i0;
  int& ix = sve.ix;
  float& xr = sve.xr;
  if (is_called_first_time) {
    i0 = 0;
  }
  //C                ********This gives the eikonal function from a table
  //C                        calculated in the first call
  //Clin-10/29/02 unsaved FR causes wrong values for ROMG with f77 compiler:
  //Cc        SAVE FR
  //C
  if (i0 != 0) {
    goto statement_100;
  }
  FEM_DO_SAFE(i, 1, 1001) {
    xr = (i - 1) * 0.01f;
    fr(i - 1) = omg0(cmn, xr);
  }
  statement_100:
  i0 = 1;
  if (x >= 10.0f) {
    return_value = 0.0f;
    return return_value;
  }
  ix = fem::fint(x * 100);
  return_value = (fr(ix) * ((ix + 1) * 0.01f - x) + fr(ix + 1) * (x -
    ix * 0.01f)) / 0.01f;
  return return_value;
}

struct hijcsc_save
{
  float bb;
  float bx;
  float by;
  float bz;
  float dis;
  float dpp1;
  float dpp2;
  float dpt1;
  float dpt2;
  float dx;
  float dy;
  float dz;
  float gs;
  float gs0;
  int i;
  int k;
  float pabs;
  arr<float> psc1;
  arr<float> psc2;
  float r2;

  hijcsc_save() :
    bb(fem::float0),
    bx(fem::float0),
    by(fem::float0),
    bz(fem::float0),
    dis(fem::float0),
    dpp1(fem::float0),
    dpp2(fem::float0),
    dpt1(fem::float0),
    dpt2(fem::float0),
    dx(fem::float0),
    dy(fem::float0),
    dz(fem::float0),
    gs(fem::float0),
    gs0(fem::float0),
    i(fem::int0),
    k(fem::int0),
    pabs(fem::float0),
    psc1(dimension(5), fem::fill0),
    psc2(dimension(5), fem::fill0),
    r2(fem::float0)
  {}
};

//C
//C*******************************************************************
//C        This subroutine performs elastic scatterings and possible
//C        elastic cascading within their own nuclei
//C*******************************************************************
void
hijcsc(
  common& cmn,
  int const& jp,
  int const& jt)
{
  FEM_CMN_SVE(hijcsc);
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  int& nseed = cmn.nseed;
  arr_ref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_ref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  //
  float& bb = sve.bb;
  float& bx = sve.bx;
  float& by = sve.by;
  float& bz = sve.bz;
  float& dis = sve.dis;
  float& dpp1 = sve.dpp1;
  float& dpp2 = sve.dpp2;
  float& dpt1 = sve.dpt1;
  float& dpt2 = sve.dpt2;
  float& dx = sve.dx;
  float& dy = sve.dy;
  float& dz = sve.dz;
  float& gs = sve.gs;
  float& gs0 = sve.gs0;
  int& i = sve.i;
  int& k = sve.k;
  float& pabs = sve.pabs;
  arr_ref<float> psc1(sve.psc1, dimension(5));
  arr_ref<float> psc2(sve.psc2, dimension(5));
  float& r2 = sve.r2;
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /RNDF77/
  //Cc      SAVE /HSTRNG/
  if (jp == 0 || jt == 0) {
    goto statement_25;
  }
  FEM_DO_SAFE(i, 1, 5) {
    psc1(i) = pp(jp, i);
    psc2(i) = pt(jt, i);
  }
  hijels(cmn, psc1, psc2);
  dpp1 = psc1(1) - pp(jp, 1);
  dpp2 = psc1(2) - pp(jp, 2);
  dpt1 = psc2(1) - pt(jt, 1);
  dpt2 = psc2(2) - pt(jt, 2);
  pp(jp, 6) += dpp1 / 2.0f;
  pp(jp, 7) += dpp2 / 2.0f;
  pp(jp, 8) += dpp1 / 2.0f;
  pp(jp, 9) += dpp2 / 2.0f;
  pt(jt, 6) += dpt1 / 2.0f;
  pt(jt, 7) += dpt2 / 2.0f;
  pt(jt, 8) += dpt1 / 2.0f;
  pt(jt, 9) += dpt2 / 2.0f;
  FEM_DO_SAFE(i, 1, 4) {
    pp(jp, i) = psc1(i);
    pt(jt, i) = psc2(i);
  }
  nfp(jp, 5) = fem::max(1, nfp(jp, 5));
  nft(jt, 5) = fem::max(1, nft(jt, 5));
  //C                ********Perform elastic scattering between JP and JT
  return;
  //C                ********The following is for possible elastic cascade
  //C
  statement_25:
  if (jp == 0) {
    goto statement_45;
  }
  pabs = fem::sqrt(fem::pow2(pp(jp, 1)) + fem::pow2(pp(jp, 2)) +
    fem::pow2(pp(jp, 3)));
  bx = pp(jp, 1) / pabs;
  by = pp(jp, 2) / pabs;
  bz = pp(jp, 3) / pabs;
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    if (i == jp) {
      goto statement_40;
    }
    dx = yp(1, i) - yp(1, jp);
    dy = yp(2, i) - yp(2, jp);
    dz = yp(3, i) - yp(3, jp);
    dis = dx * bx + dy * by + dz * bz;
    if (dis <= 0) {
      goto statement_40;
    }
    bb = fem::pow2(dx) + fem::pow2(dy) + fem::pow2(dz) - fem::pow2(dis);
    r2 = bb * hipr1(40) / hipr1(31) / 0.1f;
    //C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
    gs = 1.0f - fem::pow2(fem::exp(-(hipr1(30) + hint1(11)) / hipr1(
      31) / 2.0f * romg(cmn, r2)));
    gs0 = 1.0f - fem::pow2(fem::exp(-(hipr1(30) + hint1(11)) / hipr1(
      31) / 2.0f * romg(cmn, 0.0f)));
    if (ranart(nseed) > gs / gs0) {
      goto statement_40;
    }
    FEM_DO_SAFE(k, 1, 5) {
      psc1(k) = pp(jp, k);
      psc2(k) = pp(i, k);
    }
    hijels(cmn, psc1, psc2);
    dpp1 = psc1(1) - pp(jp, 1);
    dpp2 = psc1(2) - pp(jp, 2);
    dpt1 = psc2(1) - pp(i, 1);
    dpt2 = psc2(2) - pp(i, 2);
    pp(jp, 6) += dpp1 / 2.0f;
    pp(jp, 7) += dpp2 / 2.0f;
    pp(jp, 8) += dpp1 / 2.0f;
    pp(jp, 9) += dpp2 / 2.0f;
    pp(i, 6) += dpt1 / 2.0f;
    pp(i, 7) += dpt2 / 2.0f;
    pp(i, 8) += dpt1 / 2.0f;
    pp(i, 9) += dpt2 / 2.0f;
    FEM_DO_SAFE(k, 1, 5) {
      pp(jp, k) = psc1(k);
      pp(i, k) = psc2(k);
    }
    nfp(i, 5) = fem::max(1, nfp(i, 5));
    goto statement_45;
    statement_40:;
  }
  statement_45:
  if (jt == 0) {
    goto statement_80;
  }
  //Clin 50        PABS=SQRT(PT(JT,1)**2+PT(JT,2)**2+PT(JT,3)**2)
  pabs = fem::sqrt(fem::pow2(pt(jt, 1)) + fem::pow2(pt(jt, 2)) +
    fem::pow2(pt(jt, 3)));
  bx = pt(jt, 1) / pabs;
  by = pt(jt, 2) / pabs;
  bz = pt(jt, 3) / pabs;
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    if (i == jt) {
      goto statement_70;
    }
    dx = yt(1, i) - yt(1, jt);
    dy = yt(2, i) - yt(2, jt);
    dz = yt(3, i) - yt(3, jt);
    dis = dx * bx + dy * by + dz * bz;
    if (dis <= 0) {
      goto statement_70;
    }
    bb = fem::pow2(dx) + fem::pow2(dy) + fem::pow2(dz) - fem::pow2(dis);
    r2 = bb * hipr1(40) / hipr1(31) / 0.1f;
    //C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
    gs = fem::pow2((1.0f - fem::exp(-(hipr1(30) + hint1(11)) / hipr1(
      31) / 2.0f * romg(cmn, r2))));
    gs0 = fem::pow2((1.0f - fem::exp(-(hipr1(30) + hint1(11)) / hipr1(
      31) / 2.0f * romg(cmn, 0.0f))));
    if (ranart(nseed) > gs / gs0) {
      goto statement_70;
    }
    FEM_DO_SAFE(k, 1, 5) {
      psc1(k) = pt(jt, k);
      psc2(k) = pt(i, k);
    }
    hijels(cmn, psc1, psc2);
    dpp1 = psc1(1) - pt(jt, 1);
    dpp2 = psc1(2) - pt(jt, 2);
    dpt1 = psc2(1) - pt(i, 1);
    dpt2 = psc2(2) - pt(i, 2);
    pt(jt, 6) += dpp1 / 2.0f;
    pt(jt, 7) += dpp2 / 2.0f;
    pt(jt, 8) += dpp1 / 2.0f;
    pt(jt, 9) += dpp2 / 2.0f;
    pt(i, 6) += dpt1 / 2.0f;
    pt(i, 7) += dpt2 / 2.0f;
    pt(i, 8) += dpt1 / 2.0f;
    pt(i, 9) += dpt2 / 2.0f;
    FEM_DO_SAFE(k, 1, 5) {
      pt(jt, k) = psc1(k);
      pt(i, k) = psc2(k);
    }
    nft(i, 5) = fem::max(1, nft(i, 5));
    goto statement_80;
    statement_70:;
  }
  statement_80:;
}

struct hirnd2_save
{
  int j;
  int jl;
  int jm;
  int jmax;
  int jmin;
  int ju;
  float rx;

  hirnd2_save() :
    j(fem::int0),
    jl(fem::int0),
    jm(fem::int0),
    jmax(fem::int0),
    jmin(fem::int0),
    ju(fem::int0),
    rx(fem::float0)
  {}
};

//C
//C        This generate random number between XMIN and XMAX
float
hirnd2(
  common& cmn,
  int const& i,
  float& xmin,
  float& xmax)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(hirnd2);
  arr_cref<float, 2> rr(cmn.rr, dimension(10, 201));
  arr_cref<float, 2> xx(cmn.xx, dimension(10, 201));
  //
  int& j = sve.j;
  int& jl = sve.jl;
  int& jm = sve.jm;
  int& jmax = sve.jmax;
  int& jmin = sve.jmin;
  int& ju = sve.ju;
  float& rx = sve.rx;
  //Cc      SAVE /HIJHB/
  //Cc      SAVE /RNDF77/
  if (xmin < xx(i, 1)) {
    xmin = xx(i, 1);
  }
  if (xmax > xx(i, 201)) {
    xmax = xx(i, 201);
  }
  jmin = 1 + fem::fint(200 * (xmin - xx(i, 1)) / (xx(i, 201) - xx(i, 1)));
  jmax = 1 + fem::fint(200 * (xmax - xx(i, 1)) / (xx(i, 201) - xx(i, 1)));
  rx = rr(i, jmin) + (rr(i, jmax) - rr(i, jmin)) * ranart(cmn.nseed);
  jl = 0;
  ju = 202;
  statement_10:
  if (ju - jl > 1) {
    jm = (ju + jl) / 2;
    if ((rr(i, 201) > rr(i, 1)) == (rx > rr(i, jm))) {
      jl = jm;
    }
    else {
      ju = jm;
    }
    goto statement_10;
  }
  j = jl;
  if (j < 1) {
    j = 1;
  }
  if (j >= 201) {
    j = 200;
  }
  return_value = (xx(i, j) + xx(i, j + 1)) / 2.0f;
  return return_value;
}

struct hijsft_save
{
  float ampd;
  float ampx;
  float amq;
  float amtd;
  float amtx;
  float amx;
  float bb;
  float bb1;
  float bb2;
  float bx;
  float by;
  float cc;
  float cthep;
  float cthet;
  float d1;
  float d2;
  float dd;
  float dd1;
  float dd2;
  float dd3;
  float dd4;
  float dp1;
  float dp2;
  float dpd;
  float dpe1;
  float dpe2;
  float dpkc11;
  float dpkc12;
  float dpkc21;
  float dpkc22;
  float dpm0;
  float dpn;
  float dpx;
  float dpx1;
  float dpx2;
  float dpy1;
  float dpy2;
  float dpz1;
  float dpz2;
  float dtd;
  float dtm0;
  float dtn;
  float dtx;
  float epm;
  float epmprm;
  float epp;
  float eppprm;
  float etm;
  float etmprm;
  float etp;
  float etpprm;
  int isng;
  int jsb;
  int kcdip;
  int kcdit;
  int miss;
  int miss4;
  int nfp3;
  int nfp5;
  int nft3;
  int nft5;
  int nsb;
  float phi;
  float phi0;
  float phi1;
  float phi2;
  float pkc;
  float pkc1;
  float pkc11;
  float pkc12;
  float pkc2;
  float pkc21;
  float pkc22;
  float pkcmx;
  float ppjet;
  float psb;
  float ptjet;
  float ptp02;
  float ptt02;
  float r1;
  float r2;
  float sdd;
  float snn;
  float spdtn;
  float spdtx;
  float spntd;
  float spntx;
  float spxtd;
  float spxtn;
  float swptd;
  float swptn;
  float swptx;
  float sxx;
  float x1;
  float x2;
  float xmax;
  float xmax1;
  float xmax2;
  float xmaxhi;
  float xmin;
  float xmin1;
  float xmin2;
  float xminhi;
  float xp0;
  float xt0;
  float xx1;
  float xx2;
  float xxp;
  float xxt;
  float yp0;
  float yt0;

  hijsft_save() :
    ampd(fem::float0),
    ampx(fem::float0),
    amq(fem::float0),
    amtd(fem::float0),
    amtx(fem::float0),
    amx(fem::float0),
    bb(fem::float0),
    bb1(fem::float0),
    bb2(fem::float0),
    bx(fem::float0),
    by(fem::float0),
    cc(fem::float0),
    cthep(fem::float0),
    cthet(fem::float0),
    d1(fem::float0),
    d2(fem::float0),
    dd(fem::float0),
    dd1(fem::float0),
    dd2(fem::float0),
    dd3(fem::float0),
    dd4(fem::float0),
    dp1(fem::float0),
    dp2(fem::float0),
    dpd(fem::float0),
    dpe1(fem::float0),
    dpe2(fem::float0),
    dpkc11(fem::float0),
    dpkc12(fem::float0),
    dpkc21(fem::float0),
    dpkc22(fem::float0),
    dpm0(fem::float0),
    dpn(fem::float0),
    dpx(fem::float0),
    dpx1(fem::float0),
    dpx2(fem::float0),
    dpy1(fem::float0),
    dpy2(fem::float0),
    dpz1(fem::float0),
    dpz2(fem::float0),
    dtd(fem::float0),
    dtm0(fem::float0),
    dtn(fem::float0),
    dtx(fem::float0),
    epm(fem::float0),
    epmprm(fem::float0),
    epp(fem::float0),
    eppprm(fem::float0),
    etm(fem::float0),
    etmprm(fem::float0),
    etp(fem::float0),
    etpprm(fem::float0),
    isng(fem::int0),
    jsb(fem::int0),
    kcdip(fem::int0),
    kcdit(fem::int0),
    miss(fem::int0),
    miss4(fem::int0),
    nfp3(fem::int0),
    nfp5(fem::int0),
    nft3(fem::int0),
    nft5(fem::int0),
    nsb(fem::int0),
    phi(fem::float0),
    phi0(fem::float0),
    phi1(fem::float0),
    phi2(fem::float0),
    pkc(fem::float0),
    pkc1(fem::float0),
    pkc11(fem::float0),
    pkc12(fem::float0),
    pkc2(fem::float0),
    pkc21(fem::float0),
    pkc22(fem::float0),
    pkcmx(fem::float0),
    ppjet(fem::float0),
    psb(fem::float0),
    ptjet(fem::float0),
    ptp02(fem::float0),
    ptt02(fem::float0),
    r1(fem::float0),
    r2(fem::float0),
    sdd(fem::float0),
    snn(fem::float0),
    spdtn(fem::float0),
    spdtx(fem::float0),
    spntd(fem::float0),
    spntx(fem::float0),
    spxtd(fem::float0),
    spxtn(fem::float0),
    swptd(fem::float0),
    swptn(fem::float0),
    swptx(fem::float0),
    sxx(fem::float0),
    x1(fem::float0),
    x2(fem::float0),
    xmax(fem::float0),
    xmax1(fem::float0),
    xmax2(fem::float0),
    xmaxhi(fem::float0),
    xmin(fem::float0),
    xmin1(fem::float0),
    xmin2(fem::float0),
    xminhi(fem::float0),
    xp0(fem::float0),
    xt0(fem::float0),
    xx1(fem::float0),
    xx2(fem::float0),
    xxp(fem::float0),
    xxt(fem::float0),
    yp0(fem::float0),
    yt0(fem::float0)
  {}
};

//C
//C*******************************************************************
//C                                                                      *
//C                Subroutine HIJSFT                                   *
//C                                                                   *
//C  Scatter two excited strings, JP from proj and JT from target    *
//C*******************************************************************
void
hijsft(
  common& cmn,
  int const& jp,
  int const& jt,
  int const& jout,
  int& ierror)
{
  FEM_CMN_SVE(hijsft);
  common_write write(cmn);
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_ref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<float> hidat(cmn.hidat, dimension(10));
  int& nseed = cmn.nseed;
  arr_ref<int> npj(cmn.npj, dimension(300));
  arr_ref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_ref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_ref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_ref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_ref<float, 2> pjpe(cmn.pjpe, dimension(300, 500));
  arr_ref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_ref<int> ntj(cmn.ntj, dimension(300));
  arr_ref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_ref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_ref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_ref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_ref<float, 2> pjte(cmn.pjte, dimension(300, 500));
  arr_ref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  arr_ref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_ref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  float& ampn = cmn.ampn;
  float& amtn = cmn.amtn;
  float& amp0 = cmn.amp0;
  float& amt0 = cmn.amt0;
  int& nfdp = cmn.nfdp;
  int& nfdt = cmn.nfdt;
  float& wp = cmn.wp;
  float& wm = cmn.wm;
  float& sw = cmn.sw;
  float& dpkc1 = cmn.dpkc1;
  float& dpkc2 = cmn.dpkc2;
  float& pp11 = cmn.pp11;
  float& pp12 = cmn.pp12;
  float& pt11 = cmn.pt11;
  float& pt12 = cmn.pt12;
  float& ptp2 = cmn.ptp2;
  float& ptt2 = cmn.ptt2;
  //
  float& ampd = sve.ampd;
  float& ampx = sve.ampx;
  float& amq = sve.amq;
  float& amtd = sve.amtd;
  float& amtx = sve.amtx;
  float& amx = sve.amx;
  float& bb = sve.bb;
  float& bb1 = sve.bb1;
  float& bb2 = sve.bb2;
  float& bx = sve.bx;
  float& by = sve.by;
  float& cc = sve.cc;
  float& cthep = sve.cthep;
  float& cthet = sve.cthet;
  float& d1 = sve.d1;
  float& d2 = sve.d2;
  float& dd = sve.dd;
  float& dd1 = sve.dd1;
  float& dd2 = sve.dd2;
  float& dd3 = sve.dd3;
  float& dd4 = sve.dd4;
  float& dp1 = sve.dp1;
  float& dp2 = sve.dp2;
  float& dpd = sve.dpd;
  float& dpe1 = sve.dpe1;
  float& dpe2 = sve.dpe2;
  float& dpkc11 = sve.dpkc11;
  float& dpkc12 = sve.dpkc12;
  float& dpkc21 = sve.dpkc21;
  float& dpkc22 = sve.dpkc22;
  float& dpm0 = sve.dpm0;
  float& dpn = sve.dpn;
  float& dpx = sve.dpx;
  float& dpx1 = sve.dpx1;
  float& dpx2 = sve.dpx2;
  float& dpy1 = sve.dpy1;
  float& dpy2 = sve.dpy2;
  float& dpz1 = sve.dpz1;
  float& dpz2 = sve.dpz2;
  float& dtd = sve.dtd;
  float& dtm0 = sve.dtm0;
  float& dtn = sve.dtn;
  float& dtx = sve.dtx;
  float& epm = sve.epm;
  float& epmprm = sve.epmprm;
  float& epp = sve.epp;
  float& eppprm = sve.eppprm;
  float& etm = sve.etm;
  float& etmprm = sve.etmprm;
  float& etp = sve.etp;
  float& etpprm = sve.etpprm;
  int& isng = sve.isng;
  int& jsb = sve.jsb;
  int& kcdip = sve.kcdip;
  int& kcdit = sve.kcdit;
  int& miss = sve.miss;
  int& miss4 = sve.miss4;
  int& nfp3 = sve.nfp3;
  int& nfp5 = sve.nfp5;
  int& nft3 = sve.nft3;
  int& nft5 = sve.nft5;
  int& nsb = sve.nsb;
  float& phi = sve.phi;
  float& phi0 = sve.phi0;
  float& phi1 = sve.phi1;
  float& phi2 = sve.phi2;
  float& pkc = sve.pkc;
  float& pkc1 = sve.pkc1;
  float& pkc11 = sve.pkc11;
  float& pkc12 = sve.pkc12;
  float& pkc2 = sve.pkc2;
  float& pkc21 = sve.pkc21;
  float& pkc22 = sve.pkc22;
  float& pkcmx = sve.pkcmx;
  float& ppjet = sve.ppjet;
  float& psb = sve.psb;
  float& ptjet = sve.ptjet;
  float& ptp02 = sve.ptp02;
  float& ptt02 = sve.ptt02;
  float& r1 = sve.r1;
  float& r2 = sve.r2;
  float& sdd = sve.sdd;
  float& snn = sve.snn;
  float& spdtn = sve.spdtn;
  float& spdtx = sve.spdtx;
  float& spntd = sve.spntd;
  float& spntx = sve.spntx;
  float& spxtd = sve.spxtd;
  float& spxtn = sve.spxtn;
  float& swptd = sve.swptd;
  float& swptn = sve.swptn;
  float& swptx = sve.swptx;
  float& sxx = sve.sxx;
  float& x1 = sve.x1;
  float& x2 = sve.x2;
  float& xmax = sve.xmax;
  float& xmax1 = sve.xmax1;
  float& xmax2 = sve.xmax2;
  float& xmaxhi = sve.xmaxhi;
  float& xmin = sve.xmin;
  float& xmin1 = sve.xmin1;
  float& xmin2 = sve.xmin2;
  float& xminhi = sve.xminhi;
  float& xp0 = sve.xp0;
  float& xt0 = sve.xt0;
  float& xx1 = sve.xx1;
  float& xx2 = sve.xx2;
  float& xxp = sve.xxp;
  float& xxt = sve.xxt;
  float& yp0 = sve.yp0;
  float& yt0 = sve.yt0;
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HIJDAT/
  //Cc      SAVE /RNDF77/
  //Cc      SAVE /HJJET1/
  //Clin-4/25/01
  //C        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
  //C     &                K2SG(900,100),PXSG(900,100),PYSG(900,100),
  //C     &                PZSG(900,100),PESG(900,100),PMSG(900,100)
  //Cc      SAVE /HJJET2/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /DPMCM1/
  //Cc      SAVE /DPMCM2/
  //C*******************************************************************
  //C        JOUT-> the number
  //C        of hard scatterings preceding this soft collision.
  //C       IHNT2(13)-> 1=
  //C        double diffrac 2=single diffrac, 3=non-single diffrac.
  //C*******************************************************************
  ierror = 0;
  cmn.jjp = jp;
  cmn.jjt = jt;
  cmn.ndpm = 0;
  //C        IOPMAIN=0
  if (jp > ihnt2(1) || jt > ihnt2(3)) {
    return;
  }
  //C
  epp = pp(jp, 4) + pp(jp, 3);
  epm = pp(jp, 4) - pp(jp, 3);
  etp = pt(jt, 4) + pt(jt, 3);
  etm = pt(jt, 4) - pt(jt, 3);
  //C
  wp = epp + etp;
  wm = epm + etm;
  sw = wp * wm;
  //C                ********total W+,W- and center-of-mass energy
  //C
  if (wp < 0.0f || wm < 0.0f) {
    goto statement_1000;
  }
  //C
  if (jout == 0) {
    if (epp < 0.0f) {
      goto statement_1000;
    }
    if (epm < 0.0f) {
      goto statement_1000;
    }
    if (etp < 0.0f) {
      goto statement_1000;
    }
    if (etm < 0.0f) {
      goto statement_1000;
    }
    if (epp / (epm + 0.01f) <= etp / (etm + 0.01f)) {
      return;
    }
  }
  //C                ********For strings which does not follow a jet-prod,
  //C                        scatter only if Ycm(JP)>Ycm(JT). When jets
  //C                        are produced just before this collision
  //C                        this requirement has already be enforced
  //C                        (see SUBROUTINE HIJHRD)
  ihnt2(11) = jp;
  ihnt2(12) = jt;
  //C
  miss = 0;
  pkc1 = 0.0f;
  pkc2 = 0.0f;
  pkc11 = 0.0f;
  pkc12 = 0.0f;
  pkc21 = 0.0f;
  pkc22 = 0.0f;
  dpkc11 = 0.0f;
  dpkc12 = 0.0f;
  dpkc21 = 0.0f;
  dpkc22 = 0.0f;
  if (nfp(jp, 10) == 1 || nft(jt, 10) == 1) {
    if (nfp(jp, 10) == 1) {
      phi1 = ulangl(pp(jp, 10), pp(jp, 11));
      ppjet = fem::sqrt(fem::pow2(pp(jp, 10)) + fem::pow2(pp(jp, 11)));
      pkc1 = ppjet;
      pkc11 = pp(jp, 10);
      pkc12 = pp(jp, 11);
    }
    if (nft(jt, 10) == 1) {
      phi2 = ulangl(pt(jt, 10), pt(jt, 11));
      ptjet = fem::sqrt(fem::pow2(pt(jt, 10)) + fem::pow2(pt(jt, 11)));
      pkc2 = ptjet;
      pkc21 = pt(jt, 10);
      pkc22 = pt(jt, 11);
    }
    if (ihpr2(4) > 0 && ihnt2(1) > 1 && ihnt2(3) > 1) {
      if (nfp(jp, 10) == 0) {
        phi = -phi2;
      }
      else if (nft(jt, 10) == 0) {
        phi = phi1;
      }
      else {
        phi = (phi1 + phi2 - hipr1(40)) / 2.0f;
      }
      bx = hint1(19) * fem::cos(hint1(20));
      by = hint1(19) * fem::sin(hint1(20));
      xp0 = yp(1, jp);
      yp0 = yp(2, jp);
      xt0 = yt(1, jt) + bx;
      yt0 = yt(2, jt) + by;
      r1 = fem::max(1.2f * fem::pow(ihnt2(1), 0.3333333f), fem::sqrt(
        fem::pow2(xp0) + fem::pow2(yp0)));
      r2 = fem::max(1.2f * fem::pow(ihnt2(3), 0.3333333f), fem::sqrt(
        fem::pow2((xt0 - bx)) + fem::pow2((yt0 - by))));
      if (fem::abs(fem::cos(phi)) < 1.0e-5f) {
        dd1 = r1;
        dd2 = r1;
        dd3 = fem::abs(by + fem::sqrt(fem::pow2(r2) - fem::pow2((
          xp0 - bx))) - yp0);
        dd4 = fem::abs(by - fem::sqrt(fem::pow2(r2) - fem::pow2((
          xp0 - bx))) - yp0);
        goto statement_5;
      }
      bb = 2.0f * fem::sin(phi) * (fem::cos(phi) * yp0 - fem::sin(phi) * xp0);
      cc = (fem::pow2(yp0) - fem::pow2(r1)) * fem::pow2(fem::cos(
        phi)) + xp0 * fem::sin(phi) * (xp0 * fem::sin(phi) - 2.0f *
        yp0 * fem::cos(phi));
      dd = fem::pow2(bb) - 4.0f * cc;
      if (dd < 0.0f) {
        goto statement_10;
      }
      xx1 = (-bb + fem::sqrt(dd)) / 2.0f;
      xx2 = (-bb - fem::sqrt(dd)) / 2.0f;
      dd1 = fem::abs((xx1 - xp0) / fem::cos(phi));
      dd2 = fem::abs((xx2 - xp0) / fem::cos(phi));
      //C
      bb = 2.0f * fem::sin(phi) * (fem::cos(phi) * (yt0 - by) -
        fem::sin(phi) * xt0) - 2.0f * bx;
      cc = (fem::pow2(bx) + fem::pow2((yt0 - by)) - fem::pow2(r2)) *
        fem::pow2(fem::cos(phi)) + xt0 * fem::sin(phi) * (xt0 *
        fem::sin(phi) - 2.0f * fem::cos(phi) * (yt0 - by)) - 2.0f *
        bx * fem::sin(phi) * (fem::cos(phi) * (yt0 - by) - fem::sin(
        phi) * xt0);
      dd = fem::pow2(bb) - 4.0f * cc;
      if (dd < 0.0f) {
        goto statement_10;
      }
      xx1 = (-bb + fem::sqrt(dd)) / 2.0f;
      xx2 = (-bb - fem::sqrt(dd)) / 2.0f;
      dd3 = fem::abs((xx1 - xt0) / fem::cos(phi));
      dd4 = fem::abs((xx2 - xt0) / fem::cos(phi));
      //C
      statement_5:
      dd1 = fem::min(dd1, dd3);
      dd2 = fem::min(dd2, dd4);
      if (dd1 < hipr1(13)) {
        dd1 = 0.0f;
      }
      if (dd2 < hipr1(13)) {
        dd2 = 0.0f;
      }
      if (nfp(jp, 10) == 1 && ppjet > hipr1(11)) {
        dp1 = dd1 * hipr1(14) / 2.0f;
        dp1 = fem::min(dp1, ppjet - hipr1(11));
        pkc1 = ppjet - dp1;
        dpx1 = fem::cos(phi1) * dp1;
        dpy1 = fem::sin(phi1) * dp1;
        pkc11 = pp(jp, 10) - dpx1;
        pkc12 = pp(jp, 11) - dpy1;
        if (dp1 > 0.0f) {
          cthep = pp(jp, 12) / fem::sqrt(fem::pow2(pp(jp, 12)) +
            fem::pow2(ppjet));
          dpz1 = dp1 * cthep / fem::sqrt(1.0f - fem::pow2(cthep));
          dpe1 = fem::sqrt(fem::pow2(dpx1) + fem::pow2(dpy1) + fem::pow2(dpz1));
          eppprm = pp(jp, 4) + pp(jp, 3) - dpe1 - dpz1;
          epmprm = pp(jp, 4) - pp(jp, 3) - dpe1 + dpz1;
          if (eppprm <= 0.0f || epmprm <= 0.0f) {
            goto statement_15;
          }
          epp = eppprm;
          epm = epmprm;
          pp(jp, 10) = pkc11;
          pp(jp, 11) = pkc12;
          npj(jp)++;
          kfpj(jp, npj(jp)) = 21;
          pjpx(jp, npj(jp)) = dpx1;
          pjpy(jp, npj(jp)) = dpy1;
          pjpz(jp, npj(jp)) = dpz1;
          pjpe(jp, npj(jp)) = dpe1;
          pjpm(jp, npj(jp)) = 0.0f;
          pp(jp, 3) = pp(jp, 3) - dpz1;
          pp(jp, 4) = pp(jp, 4) - dpe1;
        }
      }
      statement_15:
      if (nft(jt, 10) == 1 && ptjet > hipr1(11)) {
        dp2 = dd2 * hipr1(14) / 2.0f;
        dp2 = fem::min(dp2, ptjet - hipr1(11));
        pkc2 = ptjet - dp2;
        dpx2 = fem::cos(phi2) * dp2;
        dpy2 = fem::sin(phi2) * dp2;
        pkc21 = pt(jt, 10) - dpx2;
        pkc22 = pt(jt, 11) - dpy2;
        if (dp2 > 0.0f) {
          cthet = pt(jt, 12) / fem::sqrt(fem::pow2(pt(jt, 12)) +
            fem::pow2(ptjet));
          dpz2 = dp2 * cthet / fem::sqrt(1.0f - fem::pow2(cthet));
          dpe2 = fem::sqrt(fem::pow2(dpx2) + fem::pow2(dpy2) + fem::pow2(dpz2));
          etpprm = pt(jt, 4) + pt(jt, 3) - dpe2 - dpz2;
          etmprm = pt(jt, 4) - pt(jt, 3) - dpe2 + dpz2;
          if (etpprm <= 0.0f || etmprm <= 0.0f) {
            goto statement_16;
          }
          etp = etpprm;
          etm = etmprm;
          pt(jt, 10) = pkc21;
          pt(jt, 11) = pkc22;
          ntj(jt)++;
          kftj(jt, ntj(jt)) = 21;
          pjtx(jt, ntj(jt)) = dpx2;
          pjty(jt, ntj(jt)) = dpy2;
          pjtz(jt, ntj(jt)) = dpz2;
          pjte(jt, ntj(jt)) = dpe2;
          pjtm(jt, ntj(jt)) = 0.0f;
          pt(jt, 3) = pt(jt, 3) - dpz2;
          pt(jt, 4) = pt(jt, 4) - dpe2;
        }
      }
      statement_16:
      dpkc11 = -(pp(jp, 10) - pkc11) / 2.0f;
      dpkc12 = -(pp(jp, 11) - pkc12) / 2.0f;
      dpkc21 = -(pt(jt, 10) - pkc21) / 2.0f;
      dpkc22 = -(pt(jt, 11) - pkc22) / 2.0f;
      wp = epp + etp;
      wm = epm + etm;
      sw = wp * wm;
    }
  }
  //C                ********If jet is quenched the pt from valence quark
  //C                        hard scattering has to reduced by d*kapa
  //C
  statement_10:
  ptp02 = fem::pow2(pp(jp, 1)) + fem::pow2(pp(jp, 2));
  ptt02 = fem::pow2(pt(jt, 1)) + fem::pow2(pt(jt, 2));
  //C
  amq = fem::max(pp(jp, 14) + pp(jp, 15), pt(jt, 14) + pt(jt, 15));
  amx = hipr1(1) + amq;
  //C                ********consider mass cut-off for strings which
  //C                        must also include quark's mass
  amp0 = amx;
  dpm0 = amx;
  nfdp = 0;
  if (nfp(jp, 5) <= 2 && nfp(jp, 3) != 0) {
    amp0 = ulmass(nfp(jp, 3));
    nfdp = nfp(jp, 3) + 2 * nfp(jp, 3) / fem::abs(nfp(jp, 3));
    dpm0 = ulmass(nfdp);
    if (dpm0 <= 0.0f) {
      nfdp = nfdp - 2 * nfdp / fem::abs(nfdp);
      dpm0 = ulmass(nfdp);
    }
  }
  amt0 = amx;
  dtm0 = amx;
  nfdt = 0;
  if (nft(jt, 5) <= 2 && nft(jt, 3) != 0) {
    amt0 = ulmass(nft(jt, 3));
    nfdt = nft(jt, 3) + 2 * nft(jt, 3) / fem::abs(nft(jt, 3));
    dtm0 = ulmass(nfdt);
    if (dtm0 <= 0.0f) {
      nfdt = nfdt - 2 * nfdt / fem::abs(nfdt);
      dtm0 = ulmass(nfdt);
    }
  }
  //C
  ampn = fem::sqrt(fem::pow2(amp0) + ptp02);
  amtn = fem::sqrt(fem::pow2(amt0) + ptt02);
  snn = fem::pow2((ampn + amtn)) + 0.001f;
  //C
  if (sw < snn + 0.001f) {
    goto statement_4000;
  }
  //C                ********Scatter only if SW>SNN
  //C*****give some PT kick to the two exited strings******************
  //Clin 20        SWPTN=4.0*(MAX(AMP0,AMT0)**2+MAX(PTP02,PTT02))
  swptn = 4.0f * (fem::pow2(fem::max(amp0, amt0)) + fem::max(ptp02, ptt02));
  swptd = 4.0f * (fem::pow2(fem::max(dpm0, dtm0)) + fem::max(ptp02, ptt02));
  swptx = 4.0f * (fem::pow2(amx) + fem::max(ptp02, ptt02));
  if (sw <= swptn) {
    pkcmx = 0.0f;
  }
  else if (sw > swptn && sw <= swptd && npj(jp) == 0 && ntj(jt) == 0) {
    pkcmx = fem::sqrt(sw / 4.0f - fem::pow2(fem::max(amp0, amt0))) -
      fem::sqrt(fem::max(ptp02, ptt02));
  }
  else if (sw > swptd && sw <= swptx && npj(jp) == 0 && ntj(jt) == 0) {
    pkcmx = fem::sqrt(sw / 4.0f - fem::pow2(fem::max(dpm0, dtm0))) -
      fem::sqrt(fem::max(ptp02, ptt02));
  }
  else if (sw > swptx) {
    pkcmx = fem::sqrt(sw / 4.0f - fem::pow2(amx)) - fem::sqrt(fem::max(ptp02,
      ptt02));
  }
  //C                ********maximun PT kick
  //C*********************************************************
  //C
  if (nfp(jp, 10) == 1 || nft(jt, 10) == 1) {
    if (pkc1 > pkcmx) {
      pkc1 = pkcmx;
      pkc11 = pkc1 * fem::cos(phi1);
      pkc12 = pkc1 * fem::sin(phi1);
      dpkc11 = -(pp(jp, 10) - pkc11) / 2.0f;
      dpkc12 = -(pp(jp, 11) - pkc12) / 2.0f;
    }
    if (pkc2 > pkcmx) {
      pkc2 = pkcmx;
      pkc21 = pkc2 * fem::cos(phi2);
      pkc22 = pkc2 * fem::sin(phi2);
      dpkc21 = -(pt(jt, 10) - pkc21) / 2.0f;
      dpkc22 = -(pt(jt, 11) - pkc22) / 2.0f;
    }
    dpkc1 = dpkc11 + dpkc21;
    dpkc2 = dpkc12 + dpkc22;
    nfp(jp, 10) = -nfp(jp, 10);
    nft(jt, 10) = -nft(jt, 10);
    goto statement_40;
  }
  //C                ********If the valence quarks had a hard-collision
  //C                        the pt kick is the pt from hard-collision.
  isng = 0;
  if (ihpr2(13) != 0 && ranart(nseed) <= hidat(4)) {
    isng = 1;
  }
  if ((nfp(jp, 5) == 3 || nft(jt, 5) == 3) || (npj(jp) != 0 || nfp(jp,
      10) != 0) || (ntj(jt) != 0 || nft(jt, 10) != 0)) {
    isng = 0;
  }
  //C
  //C               ********decite whether to have single-diffractive
  if (ihpr2(5) == 0) {
    pkc = hipr1(2) * fem::sqrt(-fem::alog(1.0f - ranart(nseed) * (
      1.0f - fem::exp(-fem::pow2(pkcmx) / fem::pow2(hipr1(2))))));
    goto statement_30;
  }
  //C
  //Clin-10/28/02 get rid of argument usage mismatch in HIRND2():
  //C        PKC=HIRND2(3,0.0,PKCMX**2)
  xminhi = 0.0f;
  xmaxhi = fem::pow2(pkcmx);
  pkc = hirnd2(cmn, 3, xminhi, xmaxhi);
  //C
  pkc = fem::sqrt(pkc);
  if (pkc > hipr1(20)) {
    pkc = hipr1(2) * fem::sqrt(-fem::alog(fem::exp(-fem::pow2(hipr1(
      20)) / fem::pow2(hipr1(2))) - ranart(nseed) * (fem::exp(-fem::pow2(
      hipr1(20)) / fem::pow2(hipr1(2))) - fem::exp(-fem::pow2(pkcmx) /
      fem::pow2(hipr1(2))))));
  }
  //C
  if (isng == 1) {
    pkc = 0.65f * fem::sqrt(-fem::alog(1.0f - ranart(nseed) * (1.0f -
      fem::exp(-fem::pow2(pkcmx) / fem::pow2(0.65f)))));
  }
  //C                        ********select PT kick
  statement_30:
  phi0 = 2.0f * hipr1(40) * ranart(nseed);
  pkc11 = pkc * fem::sin(phi0);
  pkc12 = pkc * fem::cos(phi0);
  pkc21 = -pkc11;
  pkc22 = -pkc12;
  dpkc1 = 0.0f;
  dpkc2 = 0.0f;
  statement_40:
  pp11 = pp(jp, 1) + pkc11 - dpkc1;
  pp12 = pp(jp, 2) + pkc12 - dpkc2;
  pt11 = pt(jt, 1) + pkc21 - dpkc1;
  pt12 = pt(jt, 2) + pkc22 - dpkc2;
  ptp2 = fem::pow2(pp11) + fem::pow2(pp12);
  ptt2 = fem::pow2(pt11) + fem::pow2(pt12);
  //C
  ampn = fem::sqrt(fem::pow2(amp0) + ptp2);
  amtn = fem::sqrt(fem::pow2(amt0) + ptt2);
  snn = fem::pow2((ampn + amtn)) + 0.001f;
  //C***************************************
  wp = epp + etp;
  wm = epm + etm;
  sw = wp * wm;
  //C****************************************
  if (sw < snn) {
    miss++;
    if (miss <= 100) {
      pkc = 0.0f;
      goto statement_30;
    }
    if (ihpr2(10) != 0) {
      write(6, star), "Error occured in Pt kick section of HIJSFT";
    }
    goto statement_4000;
  }
  //C******************************************************************
  ampd = fem::sqrt(fem::pow2(dpm0) + ptp2);
  amtd = fem::sqrt(fem::pow2(dtm0) + ptt2);
  //C
  ampx = fem::sqrt(fem::pow2(amx) + ptp2);
  amtx = fem::sqrt(fem::pow2(amx) + ptt2);
  //C
  dpn = fem::pow2(ampn) / sw;
  dtn = fem::pow2(amtn) / sw;
  dpd = fem::pow2(ampd) / sw;
  dtd = fem::pow2(amtd) / sw;
  dpx = fem::pow2(ampx) / sw;
  dtx = fem::pow2(amtx) / sw;
  //C
  spntd = fem::pow2((ampn + amtd));
  spntx = fem::pow2((ampn + amtx));
  //C                        ********CM energy if proj=N,targ=N*
  spdtn = fem::pow2((ampd + amtn));
  spxtn = fem::pow2((ampx + amtn));
  //C                        ********CM energy if proj=N*,targ=N
  spdtx = fem::pow2((ampd + amtx));
  spxtd = fem::pow2((ampx + amtd));
  sdd = fem::pow2((ampd + amtd));
  sxx = fem::pow2((ampx + amtx));
  //C
  //C                ********CM energy if proj=delta, targ=delta
  //C****************There are many different cases**********
  //C        IF(IHPR2(15).EQ.1) GO TO 500
  //C
  //C                ********to have DPM type soft interactions
  //C
  //Clin 45        CONTINUE
  if (sw > sxx + 0.001f) {
    if (isng == 0) {
      d1 = dpx;
      d2 = dtx;
      nfp3 = 0;
      nft3 = 0;
      goto statement_400;
    }
    else {
      //C**** 5/30/1998 this is identical to the above statement. Added to
      //C**** avoid questional branching to block.
      if ((nfp(jp, 5) == 3 && nft(jt, 5) == 3) || (npj(jp) != 0 || nfp(jp,
          10) != 0) || (ntj(jt) != 0 || nft(jt, 10) != 0)) {
        d1 = dpx;
        d2 = dtx;
        nfp3 = 0;
        nft3 = 0;
        goto statement_400;
      }
      //C                ********do not allow excited strings to have
      //C                        single-diffr
      if (ranart(nseed) > 0.5f || (nft(jt, 5) > 2 || ntj(jt) != 0 || nft(jt,
          10) != 0)) {
        d1 = dpn;
        d2 = dtx;
        nfp3 = nfp(jp, 3);
        nft3 = 0;
        goto statement_220;
      }
      else {
        d1 = dpx;
        d2 = dtn;
        nfp3 = 0;
        nft3 = nft(jt, 3);
        goto statement_240;
      }
      //C                ********have single diffractive collision
    }
  }
  else if (sw > fem::max(spdtx, spxtd) + 0.001f && sw <= sxx + 0.001f) {
    if (((npj(jp) == 0 && ntj(jt) == 0 && ranart(nseed) > 0.5f) || (
        npj(jp) == 0 && ntj(jt) != 0)) && nfp(jp, 5) <= 2) {
      d1 = dpd;
      d2 = dtx;
      nfp3 = nfdp;
      nft3 = 0;
      goto statement_220;
    }
    else if (ntj(jt) == 0 && nft(jt, 5) <= 2) {
      d1 = dpx;
      d2 = dtd;
      nfp3 = 0;
      nft3 = nfdt;
      goto statement_240;
    }
    goto statement_4000;
  }
  else if (sw > fem::min(spdtx, spxtd) + 0.001f && sw <= fem::max(spdtx,
    spxtd) + 0.001f) {
    if (spdtx <= spxtd && npj(jp) == 0 && nfp(jp, 5) <= 2) {
      d1 = dpd;
      d2 = dtx;
      nfp3 = nfdp;
      nft3 = 0;
      goto statement_220;
    }
    else if (spdtx > spxtd && ntj(jt) == 0 && nft(jt, 5) <= 2) {
      d1 = dpx;
      d2 = dtd;
      nfp3 = 0;
      nft3 = nfdt;
      goto statement_240;
    }
    //C*** 5/30/1998 added to avoid questional branching to another block
    //C*** this is identical to the statement following the next ELSE IF
    if (((npj(jp) == 0 && ntj(jt) == 0 && ranart(nseed) > 0.5f) || (
        npj(jp) == 0 && ntj(jt) != 0)) && nfp(jp, 5) <= 2) {
      d1 = dpn;
      d2 = dtx;
      nfp3 = nfp(jp, 3);
      nft3 = 0;
      goto statement_220;
    }
    else if (ntj(jt) == 0 && nft(jt, 5) <= 2) {
      d1 = dpx;
      d2 = dtn;
      nfp3 = 0;
      nft3 = nft(jt, 3);
      goto statement_240;
    }
    goto statement_4000;
  }
  else if (sw > fem::max(spntx, spxtn) + 0.001f && sw <= fem::min(spdtx,
    spxtd) + 0.001f) {
    if (((npj(jp) == 0 && ntj(jt) == 0 && ranart(nseed) > 0.5f) || (
        npj(jp) == 0 && ntj(jt) != 0)) && nfp(jp, 5) <= 2) {
      d1 = dpn;
      d2 = dtx;
      nfp3 = nfp(jp, 3);
      nft3 = 0;
      goto statement_220;
    }
    else if (ntj(jt) == 0 && nft(jt, 5) <= 2) {
      d1 = dpx;
      d2 = dtn;
      nfp3 = 0;
      nft3 = nft(jt, 3);
      goto statement_240;
    }
    goto statement_4000;
  }
  else if (sw > fem::min(spntx, spxtn) + 0.001f && sw <= fem::max(spntx,
    spxtn) + 0.001f) {
    if (spntx <= spxtn && npj(jp) == 0 && nfp(jp, 5) <= 2) {
      d1 = dpn;
      d2 = dtx;
      nfp3 = nfp(jp, 3);
      nft3 = 0;
      goto statement_220;
    }
    else if (spntx > spxtn && ntj(jt) == 0 && nft(jt, 5) <= 2) {
      d1 = dpx;
      d2 = dtn;
      nfp3 = 0;
      nft3 = nft(jt, 3);
      goto statement_240;
    }
    goto statement_4000;
  }
  else if (sw <= fem::min(spntx, spxtn) + 0.001f && (npj(jp) != 0 ||
    ntj(jt) != 0)) {
    goto statement_4000;
  }
  else if (sw <= fem::min(spntx, spxtn) + 0.001f && nfp(jp, 5) > 2 && nft(jt,
    5) > 2) {
    goto statement_4000;
  }
  else if (sw > sdd + 0.001f && sw <= fem::min(spntx, spxtn) + 0.001f) {
    d1 = dpd;
    d2 = dtd;
    nfp3 = nfdp;
    nft3 = nfdt;
    goto statement_100;
  }
  else if (sw > fem::max(spntd, spdtn) + 0.001f && sw <= sdd + 0.001f) {
    if (ranart(nseed) > 0.5f) {
      d1 = dpd;
      d2 = dtn;
      nfp3 = nfdp;
      nft3 = nft(jt, 3);
      goto statement_100;
    }
    else {
      d1 = dpn;
      d2 = dtd;
      nfp3 = nfp(jp, 3);
      nft3 = nfdt;
      goto statement_100;
    }
  }
  else if (sw > fem::min(spntd, spdtn) + 0.001f && sw <= fem::max(spntd,
    spdtn) + 0.001f) {
    if (spntd > spdtn) {
      d1 = dpd;
      d2 = dtn;
      nfp3 = nfdp;
      nft3 = nft(jt, 3);
      goto statement_100;
    }
    else {
      d1 = dpn;
      d2 = dtd;
      nfp3 = nfp(jp, 3);
      nft3 = nfdt;
      goto statement_100;
    }
  }
  else if (sw <= fem::min(spntd, spdtn) + 0.001f) {
    d1 = dpn;
    d2 = dtn;
    nfp3 = nfp(jp, 3);
    nft3 = nft(jt, 3);
    goto statement_100;
  }
  write(6, star), " Error in HIJSFT: There is no path to here";
  return;
  //C
  //C***************  elastic scattering ***************
  //C        this is like elastic, both proj and targ mass
  //C        must be fixed
  //C***************************************************
  statement_100:
  nfp5 = fem::max(2, nfp(jp, 5));
  nft5 = fem::max(2, nft(jt, 5));
  bb1 = 1.0f + d1 - d2;
  bb2 = 1.0f + d2 - d1;
  if (fem::pow2(bb1) < 4.0f * d1 || fem::pow2(bb2) < 4.0f * d2) {
    miss++;
    if (miss > 100 || pkc == 0.0f) {
      goto statement_3000;
    }
    pkc = pkc * 0.5f;
    goto statement_30;
  }
  if (ranart(nseed) < 0.5f) {
    x1 = (bb1 - fem::sqrt(fem::pow2(bb1) - 4.0f * d1)) / 2.0f;
    x2 = (bb2 - fem::sqrt(fem::pow2(bb2) - 4.0f * d2)) / 2.0f;
  }
  else {
    x1 = (bb1 + fem::sqrt(fem::pow2(bb1) - 4.0f * d1)) / 2.0f;
    x2 = (bb2 + fem::sqrt(fem::pow2(bb2) - 4.0f * d2)) / 2.0f;
  }
  ihnt2(13) = 2;
  goto statement_600;
  //C
  //C********** Single diffractive ***********************
  //C either proj or targ's mass is fixed
  //C*****************************************************
  statement_220:
  nfp5 = fem::max(2, nfp(jp, 5));
  nft5 = 3;
  if (nfp3 == 0) {
    nfp5 = 3;
  }
  bb2 = 1.0f + d2 - d1;
  if (fem::pow2(bb2) < 4.0f * d2) {
    miss++;
    if (miss > 100 || pkc == 0.0f) {
      goto statement_3000;
    }
    pkc = pkc * 0.5f;
    goto statement_30;
  }
  xmin = (bb2 - fem::sqrt(fem::pow2(bb2) - 4.0f * d2)) / 2.0f;
  xmax = (bb2 + fem::sqrt(fem::pow2(bb2) - 4.0f * d2)) / 2.0f;
  miss4 = 0;
  statement_222:
  x2 = hirnd2(cmn, 6, xmin, xmax);
  x1 = d1 / (1.0f - x2);
  if (x2 * (1.0f - x1) < (d2 + 1.e-4f / sw)) {
    miss4++;
    if (miss4 <= 1000) {
      goto statement_222;
    }
    goto statement_5000;
  }
  ihnt2(13) = 2;
  goto statement_600;
  //C                        ********Fix proj mass*********
  statement_240:
  nfp5 = 3;
  nft5 = fem::max(2, nft(jt, 5));
  if (nft3 == 0) {
    nft5 = 3;
  }
  bb1 = 1.0f + d1 - d2;
  if (fem::pow2(bb1) < 4.0f * d1) {
    miss++;
    if (miss > 100 || pkc == 0.0f) {
      goto statement_3000;
    }
    pkc = pkc * 0.5f;
    goto statement_30;
  }
  xmin = (bb1 - fem::sqrt(fem::pow2(bb1) - 4.0f * d1)) / 2.0f;
  xmax = (bb1 + fem::sqrt(fem::pow2(bb1) - 4.0f * d1)) / 2.0f;
  miss4 = 0;
  statement_242:
  x1 = hirnd2(cmn, 6, xmin, xmax);
  x2 = d2 / (1.0f - x1);
  if (x1 * (1.0f - x2) < (d1 + 1.e-4f / sw)) {
    miss4++;
    if (miss4 <= 1000) {
      goto statement_242;
    }
    goto statement_5000;
  }
  ihnt2(13) = 2;
  goto statement_600;
  //C                        ********Fix targ mass*********
  //C
  //C*************non-single diffractive**********************
  //C        both proj and targ may not be fixed in mass
  //C*********************************************************
  //C
  statement_400:
  nfp5 = 3;
  nft5 = 3;
  bb1 = 1.0f + d1 - d2;
  bb2 = 1.0f + d2 - d1;
  if (fem::pow2(bb1) < 4.0f * d1 || fem::pow2(bb2) < 4.0f * d2) {
    miss++;
    if (miss > 100 || pkc == 0.0f) {
      goto statement_3000;
    }
    pkc = pkc * 0.5f;
    goto statement_30;
  }
  xmin1 = (bb1 - fem::sqrt(fem::pow2(bb1) - 4.0f * d1)) / 2.0f;
  xmax1 = (bb1 + fem::sqrt(fem::pow2(bb1) - 4.0f * d1)) / 2.0f;
  xmin2 = (bb2 - fem::sqrt(fem::pow2(bb2) - 4.0f * d2)) / 2.0f;
  xmax2 = (bb2 + fem::sqrt(fem::pow2(bb2) - 4.0f * d2)) / 2.0f;
  miss4 = 0;
  statement_410:
  x1 = hirnd2(cmn, 4, xmin1, xmax1);
  x2 = hirnd2(cmn, 4, xmin2, xmax2);
  if (nfp(jp, 5) == 3 || nft(jt, 5) == 3) {
    x1 = hirnd2(cmn, 6, xmin1, xmax1);
    x2 = hirnd2(cmn, 6, xmin2, xmax2);
  }
  //C                        ********
  if (fem::abs(nfp(jp, 1) * nfp(jp, 2)) > 1000000 || fem::abs(nfp(jp,
      1) * nfp(jp, 2)) < 100) {
    x1 = hirnd2(cmn, 5, xmin1, xmax1);
  }
  if (fem::abs(nft(jt, 1) * nft(jt, 2)) > 1000000 || fem::abs(nft(jt,
      1) * nft(jt, 2)) < 100) {
    x2 = hirnd2(cmn, 5, xmin2, xmax2);
  }
  //C        IF(IOPMAIN.EQ.3) X1=HIRND2(6,XMIN1,XMAX1)
  //C        IF(IOPMAIN.EQ.2) X2=HIRND2(6,XMIN2,XMAX2)
  //C        ********For q-qbar or (qq)-(qq)bar system use symetric
  //C                distribution, for q-(qq) or qbar-(qq)bar use
  //C                unsymetrical distribution
  //C
  if (fem::abs(nfp(jp, 1) * nfp(jp, 2)) > 1000000) {
    x1 = 1.0f - x1;
  }
  xxp = x1 * (1.0f - x2);
  xxt = x2 * (1.0f - x1);
  if (xxp < (d1 + 1.e-4f / sw) || xxt < (d2 + 1.e-4f / sw)) {
    miss4++;
    if (miss4 <= 1000) {
      goto statement_410;
    }
    goto statement_5000;
  }
  ihnt2(13) = 3;
  //C***************************************************
  statement_600:
  if (x1 * (1.0f - x2) < (fem::pow2(ampn) - 1.e-4f) / sw || x2 * (
      1.0f - x1) < (fem::pow2(amtn) - 1.e-4f) / sw) {
    miss++;
    if (miss > 100 || pkc == 0.0f) {
      goto statement_2000;
    }
    pkc = 0.0f;
    goto statement_30;
  }
  //C
  epp = (1.0f - x2) * wp;
  epm = x1 * wm;
  etp = x2 * wp;
  etm = (1.0f - x1) * wm;
  pp(jp, 3) = (epp - epm) / 2.0f;
  pp(jp, 4) = (epp + epm) / 2.0f;
  if (epp * epm - ptp2 < 0.0f) {
    goto statement_6000;
  }
  pp(jp, 5) = fem::sqrt(epp * epm - ptp2);
  nfp(jp, 3) = nfp3;
  nfp(jp, 5) = nfp5;
  //C
  pt(jt, 3) = (etp - etm) / 2.0f;
  pt(jt, 4) = (etp + etm) / 2.0f;
  if (etp * etm - ptt2 < 0.0f) {
    goto statement_6000;
  }
  pt(jt, 5) = fem::sqrt(etp * etm - ptt2);
  nft(jt, 3) = nft3;
  nft(jt, 5) = nft5;
  //C*****recoil PT from hard-inter is shared by two end-partons
  //C       so that pt=p1+p2
  pp(jp, 1) = pp11 - pkc11;
  pp(jp, 2) = pp12 - pkc12;
  //C
  kcdip = 1;
  kcdit = 1;
  if (fem::abs(nfp(jp, 1) * nfp(jp, 2)) > 1000000 || fem::abs(nfp(jp,
      1) * nfp(jp, 2)) < 100) {
    kcdip = 0;
  }
  if (fem::abs(nft(jt, 1) * nft(jt, 2)) > 1000000 || fem::abs(nft(jt,
      1) * nft(jt, 2)) < 100) {
    kcdit = 0;
  }
  if ((kcdip == 0 && ranart(nseed) < 0.5f) || (kcdip != 0 && ranart(
      nseed) < 0.5f / (1.0f + (fem::pow2(pkc11) + fem::pow2(pkc12)) /
      fem::pow2(hipr1(22))))) {
    pp(jp, 6) += (pp(jp, 1) - pp(jp, 6) - pp(jp, 8) - dpkc1) / 2.0f;
    pp(jp, 7) += (pp(jp, 2) - pp(jp, 7) - pp(jp, 9) - dpkc2) / 2.0f;
    pp(jp, 8) = (pp(jp, 1) - pp(jp, 6) - pp(jp, 8) - dpkc1) / 2.0f + pp(jp,
      8) + pkc11;
    pp(jp, 9) = (pp(jp, 2) - pp(jp, 7) - pp(jp, 9) - dpkc2) / 2.0f + pp(jp,
      9) + pkc12;
  }
  else {
    pp(jp, 8) += (pp(jp, 1) - pp(jp, 6) - pp(jp, 8) - dpkc1) / 2.0f;
    pp(jp, 9) += (pp(jp, 2) - pp(jp, 7) - pp(jp, 9) - dpkc2) / 2.0f;
    pp(jp, 6) = (pp(jp, 1) - pp(jp, 6) - pp(jp, 8) - dpkc1) / 2.0f + pp(jp,
      6) + pkc11;
    pp(jp, 7) = (pp(jp, 2) - pp(jp, 7) - pp(jp, 9) - dpkc2) / 2.0f + pp(jp,
      7) + pkc12;
  }
  pp(jp, 1) = pp(jp, 6) + pp(jp, 8);
  pp(jp, 2) = pp(jp, 7) + pp(jp, 9);
  //C                                ********pt kick for proj
  pt(jt, 1) = pt11 - pkc21;
  pt(jt, 2) = pt12 - pkc22;
  if ((kcdit == 0 && ranart(nseed) < 0.5f) || (kcdit != 0 && ranart(
      nseed) < 0.5f / (1.0f + (fem::pow2(pkc21) + fem::pow2(pkc22)) /
      fem::pow2(hipr1(22))))) {
    pt(jt, 6) += (pt(jt, 1) - pt(jt, 6) - pt(jt, 8) - dpkc1) / 2.0f;
    pt(jt, 7) += (pt(jt, 2) - pt(jt, 7) - pt(jt, 9) - dpkc2) / 2.0f;
    pt(jt, 8) = (pt(jt, 1) - pt(jt, 6) - pt(jt, 8) - dpkc1) / 2.0f + pt(jt,
      8) + pkc21;
    pt(jt, 9) = (pt(jt, 2) - pt(jt, 7) - pt(jt, 9) - dpkc2) / 2.0f + pt(jt,
      9) + pkc22;
  }
  else {
    pt(jt, 8) += (pt(jt, 1) - pt(jt, 6) - pt(jt, 8) - dpkc1) / 2.0f;
    pt(jt, 9) += (pt(jt, 2) - pt(jt, 7) - pt(jt, 9) - dpkc2) / 2.0f;
    pt(jt, 6) = (pt(jt, 1) - pt(jt, 6) - pt(jt, 8) - dpkc1) / 2.0f + pt(jt,
      6) + pkc21;
    pt(jt, 7) = (pt(jt, 2) - pt(jt, 7) - pt(jt, 9) - dpkc2) / 2.0f + pt(jt,
      7) + pkc22;
  }
  pt(jt, 1) = pt(jt, 6) + pt(jt, 8);
  pt(jt, 2) = pt(jt, 7) + pt(jt, 9);
  //C                        ********pt kick for targ
  //C
  if (npj(jp) != 0) {
    nfp(jp, 5) = 3;
  }
  if (ntj(jt) != 0) {
    nft(jt, 5) = 3;
  }
  //C                        ********jets must be connected to string
  if (epp / (epm + 0.0001f) < etp / (etm + 0.0001f) && fem::abs(nfp(jp,
      1) * nfp(jp, 2)) < 1000000) {
    FEM_DO_SAFE(jsb, 1, 15) {
      psb = pp(jp, jsb);
      pp(jp, jsb) = pt(jt, jsb);
      pt(jt, jsb) = psb;
      nsb = nfp(jp, jsb);
      nfp(jp, jsb) = nft(jt, jsb);
      nft(jt, jsb) = nsb;
    }
    //C                ********when Ycm(JP)<Ycm(JT) after the collision
    //C                        exchange the positions of the two
  }
  //C
  return;
  //C**************************************************
  statement_1000:
  ierror = 1;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), "     Fatal HIJSFT start error,abandon this event";
  write(6, star), "     PROJ E+,E-,W+", epp, epm, wp;
  write(6, star), "     TARG E+,E-,W-", etp, etm, wm;
  write(6, star), "     W+*W-, (APN+ATN)^2", sw, snn;
  return;
  statement_2000:
  ierror = 0;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), "     (2)energy partition fail,";
  write(6, star), "     HIJSFT not performed, but continue";
  write(6, star), "     MP1,MPN", x1 * (1.0f - x2) * sw, fem::pow2(ampn);
  write(6, star), "     MT2,MTN", x2 * (1.0f - x1) * sw, fem::pow2(amtn);
  return;
  statement_3000:
  ierror = 0;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), "     (3)something is wrong with the pt kick, ";
  write(6, star), "     HIJSFT not performed, but continue";
  write(6, star), "     D1=", d1, " D2=", d2, " SW=", sw;
  write(6, star), "     HISTORY NFP5=", nfp(jp, 5), " NFT5=", nft(jt, 5);
  write(6, star), "     THIS COLLISON NFP5=", nfp5, " NFT5=", nft5;
  write(6, star), "     # OF JET IN PROJ", npj(jp), " IN TARG", ntj(jt);
  return;
  statement_4000:
  ierror = 0;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), "     (4)unable to choose process, but not harmful";
  write(6, star), "     HIJSFT not performed, but continue";
  write(6, star), "     PTP=", fem::sqrt(ptp2), " PTT=", fem::sqrt(ptt2),
    " SW=", sw;
  write(6, star), "     AMCUT=", amx, " JP=", jp, " JT=", jt;
  write(6, star), "     HISTORY NFP5=", nfp(jp, 5), " NFT5=", nft(jt, 5);
  return;
  statement_5000:
  ierror = 0;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), "     energy partition failed(5),for limited try";
  write(6, star), "     HIJSFT not performed, but continue";
  write(6, star), "     NFP5=", nfp5, " NFT5=", nft5;
  write(6, star), "     D1", d1, " X1(1-X2)", x1 * (1.0f - x2);
  write(6, star), "     D2", d2, " X2(1-X1)", x2 * (1.0f - x1);
  return;
  statement_6000:
  pkc = 0.0f;
  miss++;
  if (miss < 100) {
    goto statement_30;
  }
  ierror = 1;
  if (ihpr2(10) == 0) {
    return;
  }
  write(6, star), " ERROR OCCURED, HIJSFT NOT PERFORMED";
  write(6, star), " Abort this event";
  write(6, star), "MTP,PTP2", epp * epm, ptp2, "  MTT,PTT2", etp * etm, ptt2;
}

struct hirnd_save
{
  int j;
  int jl;
  int jm;
  int ju;
  float rx;

  hirnd_save() :
    j(fem::int0),
    jl(fem::int0),
    jm(fem::int0),
    ju(fem::int0),
    rx(fem::float0)
  {}
};

float
hirnd(
  common& cmn,
  int const& i)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(hirnd);
  arr_cref<float, 2> rr(cmn.rr, dimension(10, 201));
  arr_cref<float, 2> xx(cmn.xx, dimension(10, 201));
  //
  int& j = sve.j;
  int& jl = sve.jl;
  int& jm = sve.jm;
  int& ju = sve.ju;
  float& rx = sve.rx;
  //Cc      SAVE /HIJHB/
  //Cc      SAVE /RNDF77/
  rx = ranart(cmn.nseed);
  jl = 0;
  ju = 202;
  statement_10:
  if (ju - jl > 1) {
    jm = (ju + jl) / 2;
    if ((rr(i, 201) > rr(i, 1)) == (rx > rr(i, jm))) {
      jl = jm;
    }
    else {
      ju = jm;
    }
    goto statement_10;
  }
  j = jl;
  if (j < 1) {
    j = 1;
  }
  if (j >= 201) {
    j = 200;
  }
  return_value = (xx(i, j) + xx(i, j + 1)) / 2.0f;
  return return_value;
}

struct hijing_save
{
  float aphx1;
  float aphx2;
  float b2;
  float bb;
  float bbx;
  float bby;
  float bmax;
  float bmin;
  float cx;
  float dengy;
  float dnbp;
  float dnbp1;
  float dnbp2;
  float dnbp3;
  float dnbt;
  float dnbt1;
  float dnbt2;
  float dnbt3;
  float gs;
  float gstot;
  float gstot0;
  int i;
  int i05;
  int idstr;
  int ierror;
  int ii;
  int iityp;
  int ijet;
  arr<int> ipcol;
  int isg;
  int istr;
  arr<int> itcol;
  int itest;
  int j;
  int jflg;
  int jjtp;
  int jout;
  int jp;
  int jphard;
  int jpmini;
  int jt;
  int jthard;
  int jtmini;
  arr<int> jtp;
  int kp;
  int kp2;
  int kt;
  int kt2;
  int lsg;
  int ncolt;
  int nftp;
  int nhard;
  int njet;
  int nlop;
  int nmini;
  int nmom;
  int npar;
  int npart;
  int nsbst;
  int nsbstr;
  int nstrg;
  int ntp;
  float phi;
  float r;
  float r2;
  float rantot;
  float rnd1;
  float rnd2;
  float rnd3;
  arr<float, 2> rnip;
  float rrb1;
  float rrb2;
  arr<float, 2> scip;
  arr<float, 2> sjip;
  float sx;
  float tt;
  float ttrig;
  float tts;
  float x;
  float xr;
  float xr1;
  float y1;
  float y2;
  float y3;

  hijing_save() :
    aphx1(fem::float0),
    aphx2(fem::float0),
    b2(fem::float0),
    bb(fem::float0),
    bbx(fem::float0),
    bby(fem::float0),
    bmax(fem::float0),
    bmin(fem::float0),
    cx(fem::float0),
    dengy(fem::float0),
    dnbp(fem::float0),
    dnbp1(fem::float0),
    dnbp2(fem::float0),
    dnbp3(fem::float0),
    dnbt(fem::float0),
    dnbt1(fem::float0),
    dnbt2(fem::float0),
    dnbt3(fem::float0),
    gs(fem::float0),
    gstot(fem::float0),
    gstot0(fem::float0),
    i(fem::int0),
    i05(fem::int0),
    idstr(fem::int0),
    ierror(fem::int0),
    ii(fem::int0),
    iityp(fem::int0),
    ijet(fem::int0),
    ipcol(dimension(90000), fem::fill0),
    isg(fem::int0),
    istr(fem::int0),
    itcol(dimension(90000), fem::fill0),
    itest(fem::int0),
    j(fem::int0),
    jflg(fem::int0),
    jjtp(fem::int0),
    jout(fem::int0),
    jp(fem::int0),
    jphard(fem::int0),
    jpmini(fem::int0),
    jt(fem::int0),
    jthard(fem::int0),
    jtmini(fem::int0),
    jtp(dimension(3), fem::fill0),
    kp(fem::int0),
    kp2(fem::int0),
    kt(fem::int0),
    kt2(fem::int0),
    lsg(fem::int0),
    ncolt(fem::int0),
    nftp(fem::int0),
    nhard(fem::int0),
    njet(fem::int0),
    nlop(fem::int0),
    nmini(fem::int0),
    nmom(fem::int0),
    npar(fem::int0),
    npart(fem::int0),
    nsbst(fem::int0),
    nsbstr(fem::int0),
    nstrg(fem::int0),
    ntp(fem::int0),
    phi(fem::float0),
    r(fem::float0),
    r2(fem::float0),
    rantot(fem::float0),
    rnd1(fem::float0),
    rnd2(fem::float0),
    rnd3(fem::float0),
    rnip(dimension(300, 300), fem::fill0),
    rrb1(fem::float0),
    rrb2(fem::float0),
    scip(dimension(300, 300), fem::fill0),
    sjip(dimension(300, 300), fem::fill0),
    sx(fem::float0),
    tt(fem::float0),
    ttrig(fem::float0),
    tts(fem::float0),
    x(fem::float0),
    xr(fem::float0),
    xr1(fem::float0),
    y1(fem::float0),
    y2(fem::float0),
    y3(fem::float0)
  {}
};

//C.................... hijing1.383_ampt.f
//C     Version 1.383
//C     The variables isng in HIJSFT and JL in ATTRAD were not initialized.
//C     The version initialize them. (as found by Fernando Marroquim)
//C
//C     Version 1.382
//C     Nuclear distribution for deuteron is taken as the Hulthen wave
//C     function as provided by Brian Cole (Columbia)
//Clin     used my own implementation of impact parameter
//Clin     & proton-neutron distance within a deuteron.
//C
//C     Version 1.381
//C
//C     The parameters for Wood-Saxon distribution for deuteron are
//C     constrained to give the right rms ratius 2.116 fm
//C     (R=0.0, D=0.5882)
//C
//C     Version 1.38
//C
//C     The following common block is added to record the number of elastic
//C     (NELT, NELP) and inelastic (NINT, NINP) participants
//C
//C        COMMON/HJGLBR/NELT,NINT,NELP,NINP
//C        SAVE /HJGLBR/
//C
//C     Version 1.37
//C
//C     A bug in the quenching subroutine is corrected. When calculating the
//C     distance between two wounded nucleons, the displacement of the
//C     impact parameter was not inculded. This bug was discovered by
//C     Dr. V.Uzhinskii JINR, Dubna, Russia
//C
//C     Version 1.36
//C
//C     Modification Oct. 8, 1998. In hijing, log(ran(nseed)) occasionally
//C     causes overfloat. It is modified to log(max(ran(nseed),1.0e-20)).
//C
//C     Nothing important has been changed here. A few 'garbage' has been
//C     cleaned up here, like common block HJJET3 for the sea quark strings
//C     which were originally created to implement the DPM scheme which
//C     later was abadoned in the final version. The lines which operate
//C     on these data are also deleted in the program.
//C
//C     Version 1.35
//C     There are some changes in the program: subroutine HARDJET is now
//C     consolidated with HIJHRD. HARDJET is used to re-initiate PYTHIA
//C     for the triggered hard processes. Now that is done  altogether
//C     with other normal hard processes in modified JETINI. In the new
//C     version one calls JETINI every time one calls HIJHRD. In the new
//C     version the effect of the isospin of the nucleon on hard processes,
//C     especially direct photons is correctly considered.
//C     For A+A collisions, one has to initilize pythia
//C     separately for each type of collisions, pp, pn,np and nn,
//C     or hp and hn for hA collisions. In JETINI we use the following
//C     catalogue for different types of collisions:
//C     h+h: h+h (itype=1)
//C     h+A: h+p (itype=1), h+n (itype=2)
//C     A+h: p+h (itype=1), n+h (itype=2)
//C     A+A: p+p (itype=1), p+n (itype=2), n+p (itype=3), n+n (itype=4)
//C*****************************************************************
//C
//C     Version 1.34
//C     Last modification on January 5, 1998. Two mistakes are corrected in
//C     function G. A Mistake in the subroutine Parton is also corrected.
//C     (These are pointed out by Ysushi Nara).
//C
//C       Last modifcation on April 10, 1996. To conduct final
//C       state radiation, PYTHIA reorganize the two scattered
//C       partons and their final momenta will be a little
//C       different. The summed total momenta of the partons
//C       from the final state radiation are stored in HINT1(26-29)
//C       and HINT1(36-39) which are little different from
//C       HINT1(21-24) and HINT1(41-44).
//C
//C       Version 1.33
//C
//C       Last modfication  on September 11, 1995. When HIJING and
//C       PYTHIA are initialized, the shadowing is evaluated at
//C       b=0 which is the maximum. This will cause overestimate
//C       of shadowing for peripheral interactions. To correct this
//C       problem, shadowing is set to zero when initializing. Then
//C       use these maximum  cross section without shadowing as a
//C       normalization of the Monte Carlo. This however increase
//C       the computing time. IHNT2(16) is used to indicate whether
//C       the sturcture function is called for (IHNT2(16)=1) initialization
//C       or for (IHNT2(16)=0)normal collisions simulation
//C
//C       Last modification on Aagust 28, 1994. Two bugs associate
//C       with the impact parameter dependence of the shadowing is
//C       corrected.
//C
//C       Last modification on October 14, 1994. One bug is corrected
//C       in the direct photon production option in subroutine
//C       HIJHRD.( this problem was reported by Jim Carroll and Mike Beddo).
//C       Another bug associated with keeping the decay history
//C       in the particle information is also corrected.(this problem
//C       was reported by Matt Bloomer)
//C
//C       Last modification on July 15, 1994. The option to trig on
//C       heavy quark production (charm IHPR2(18)=0 or beauty IHPR2(18)=1)
//C       is added. To do this, set IHPR2(3)=3. For inclusive production,
//C       one should reset HIPR1(10)=0.0. One can also trig larger pt
//C       QQbar production by giving HIPR1(10) a nonvanishing value.
//C       The mass of the heavy quark in the calculation of the cross
//C       section (HINT1(59)--HINT1(65)) is given by HIPR1(7) (the
//C       default is the charm mass D=1.5). We also include a separate
//C       K-factor for heavy quark and direct photon production by
//C       HIPR1(23)(D=2.0).
//C
//C       Last modification on May 24, 1994.  The option to
//C       retain the information of all particles including those
//C       who have decayed is IHPR(21)=1 (default=0). KATT(I,3) is
//C       added to contain the line number of the parent particle
//C       of the current line which is produced via a decay.
//C       KATT(I,4) is the status number of the particle: 11=particle
//C       which has decayed; 1=finally produced particle.
//C
//C       Last modification on May 24, 1994( in HIJSFT when valence quark
//C       is quenched, the following error is corrected. 1.2*IHNT2(1) -->
//C       1.2*IHNT2(1)**0.333333, 1.2*IHNT2(3) -->1.2*IHNT(3)**0.333333)
//C
//C       Last modification on March 16, 1994 (heavy flavor production
//C       processes MSUB(81)=1 MSUB(82)=1 have been switched on,
//C       charm production is the default, B-quark option is
//C       IHPR2(18), when it is switched on, charm quark is
//C       automatically off)
//C
//C       Last modification on March 23, 1994 (an error is corrected
//C       in the impact parameter dependence of the jet cross section)
//C
//C       Last modification Oct. 1993 to comply with non-vax
//C       machines' compiler
//C
//C*********************************************
//C        LAST MODIFICATION April 5, 1991
//CQUARK DISTRIBUTIOIN (1-X)**A/(X**2+C**2/S)**B
//C(A=HIPR1(44),B=HIPR1(46),C=HIPR1(45))
//C STRING FLIP, VENUS OPTION IHPR2(15)=1,IN WHICH ONE CAN HAVE ONE AND
//C TWO COLOR CHANGES, (1-W)**2,W*(1-W),W*(1-W),AND W*2, W=HIPR1(18),
//C AMONG PT DISTRIBUTION OF SEA QUARKS IS CONTROLLED BY HIPR1(42)
//C
//C        gluon jets can form a single string system
//C
//C        initial state radiation is included
//C
//C        all QCD subprocesses are included
//C
//C        direct particles production is included(currently only direct
//C                photon)
//C
//C        Effect of high P_T trigger bias on multiple jets distribution
//C
//C******************************************************************
//C                                HIJING.10                         *
//C                  Heavy Ion Jet INteraction Generator             *
//C                                   by                             *
//C                   X. N. Wang      and   M. Gyulassy              *
//C                      Lawrence Berkeley Laboratory                *
//C                                                                  *
//C******************************************************************
//C
//C******************************************************************
//C NFP(K,1),NFP(K,2)=flavor of q and di-q, NFP(K,3)=present ID of  *
//C proj, NFP(K,4) original ID of proj.  NFP(K,5)=colli status(0=no,*
//C 1=elastic,2=the diffrac one in single-diffrac,3= excited string.*
//C |NFP(K,6)| is the total # of jet production, if NFP(K,6)<0 it   *
//C can not produce jet anymore. NFP(K,10)=valence quarks scattering*
//C (0=has not been,1=is going to be, -1=has already been scattered *
//C NFP(k,11) total number of interactions this proj has suffered   *
//C PP(K,1)=PX,PP(K,2)=PY,PP(K,3)=PZ,PP(K,4)=E,PP(K,5)=M(invariant  *
//C mass), PP(K,6,7),PP(K,8,9)=transverse momentum of quark and     *
//C diquark,PP(K,10)=PT of the hard scattering between the valence  *
//C quarks; PP(K,14,15)=the mass of quark,diquark.                   *
//C******************************************************************
//C
//C****************************************************************
//C
//C        SUBROUTINE HIJING
//C
//C****************************************************************
void
hijing(
  common& cmn,
  str_cref frame,
  float const& bmin0,
  float const& bmax0)
{
  FEM_CMN_SVE(hijing);
  common_write write(cmn);
  arr_ref<float> hipr1(cmn.hipr1, dimension(100));
  arr_ref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  arr_ref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_ref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_ref<float, 2> yt(cmn.yt, dimension(3, 300));
  int& nelt = cmn.nelt;
  int& ninthj = cmn.ninthj;
  int& nelp = cmn.nelp;
  int& ninp = cmn.ninp;
  float& eatt = cmn.eatt;
  int& jatt = cmn.jatt;
  int& natt = cmn.natt;
  int& nt = cmn.nt;
  int& np = cmn.np;
  int& n0 = cmn.n0;
  int& n01 = cmn.n01;
  int& n10 = cmn.n10;
  int& n11 = cmn.n11;
  const int maxstr = 150001;
  arr_ref<int, 2> katt(cmn.katt, dimension(maxstr, 4));
  arr_ref<float, 2> patt(cmn.patt, dimension(maxstr, 4));
  arr_ref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_ref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  arr_ref<int> npj(cmn.npj, dimension(300));
  arr_ref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_ref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_ref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_ref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_ref<float, 2> pjpe(cmn.pjpe, dimension(300, 500));
  arr_ref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_ref<int> ntj(cmn.ntj, dimension(300));
  arr_ref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_ref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_ref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_ref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_ref<float, 2> pjte(cmn.pjte, dimension(300, 500));
  arr_ref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  int& nsg = cmn.nsg;
  arr_ref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_ref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_ref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_ref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_ref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_ref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_ref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_ref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  int& ndr = cmn.ndr;
  arr_cref<int> kfdr(cmn.kfdr, dimension(maxstr));
  arr_cref<float, 2> pdr(cmn.pdr, dimension(maxstr, 5));
  arr_cref<float, 2> rtdr(cmn.rtdr, dimension(maxstr, 2));
  int& nseed = cmn.nseed;
  int& n = static_cast<common_lujets&>(cmn).n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_cref<float, 2> p(cmn.p, dimension(9000, 5));
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_ref<int> itypar(cmn.itypar, dimension(maxstr));
  arr_ref<float> gxar(cmn.gxar, dimension(maxstr));
  arr_ref<float> gyar(cmn.gyar, dimension(maxstr));
  arr_ref<float> gzar(cmn.gzar, dimension(maxstr));
  arr_ref<float> ftar(cmn.ftar, dimension(maxstr));
  arr_ref<float> pxar(cmn.pxar, dimension(maxstr));
  arr_ref<float> pyar(cmn.pyar, dimension(maxstr));
  arr_ref<float> pzar(cmn.pzar, dimension(maxstr));
  arr_ref<float> pear(cmn.pear, dimension(maxstr));
  arr_ref<float> xmar(cmn.xmar, dimension(maxstr));
  int& mul = cmn.mul;
  const int maxptn = 400001;
  arr_ref<double> gx0(cmn.gx0, dimension(maxptn));
  arr_ref<double> gy0(cmn.gy0, dimension(maxptn));
  arr_ref<double> gz0(cmn.gz0, dimension(maxptn));
  arr_ref<double> ft0(cmn.ft0, dimension(maxptn));
  arr_ref<double> px0(cmn.px0, dimension(maxptn));
  arr_ref<double> py0(cmn.py0, dimension(maxptn));
  arr_ref<double> pz0(cmn.pz0, dimension(maxptn));
  arr_ref<double> e0(cmn.e0, dimension(maxptn));
  arr_ref<double> xmass0(cmn.xmass0, dimension(maxptn));
  arr_ref<int> ityp0(cmn.ityp0, dimension(maxptn));
  arr_cref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_cref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_cref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_cref<double> px5(cmn.px5, dimension(maxptn));
  arr_cref<double> py5(cmn.py5, dimension(maxptn));
  arr_cref<double> pz5(cmn.pz5, dimension(maxptn));
  arr_ref<double> e5(cmn.e5, dimension(maxptn));
  arr_cref<double> xmass5(cmn.xmass5, dimension(maxptn));
  arr_cref<int> ityp5(cmn.ityp5, dimension(maxptn));
  arr_ref<int> lstrg0(cmn.lstrg0, dimension(maxptn));
  arr_ref<int> lpart0(cmn.lpart0, dimension(maxptn));
  arr_cref<int> lstrg1(cmn.lstrg1, dimension(maxptn));
  arr_cref<int> lpart1(cmn.lpart1, dimension(maxptn));
  int& nsp = cmn.nsp;
  int& nst = cmn.nst;
  int& nsi = cmn.nsi;
  arr_cref<double> ataui(cmn.ataui, dimension(maxstr));
  arr_cref<double> zt1(cmn.zt1, dimension(maxstr));
  arr_cref<double> zt2(cmn.zt2, dimension(maxstr));
  arr_cref<double> zt3(cmn.zt3, dimension(maxstr));
  int& isoft = cmn.isoft;
  int& isflag = cmn.isflag;
  arr_ref<double, 2> pxsgs(cmn.pxsgs, dimension(maxstr, 3));
  arr_ref<double, 2> pysgs(cmn.pysgs, dimension(maxstr, 3));
  arr_ref<double, 2> pzsgs(cmn.pzsgs, dimension(maxstr, 3));
  arr_ref<double, 2> pesgs(cmn.pesgs, dimension(maxstr, 3));
  arr_ref<double, 2> pmsgs(cmn.pmsgs, dimension(maxstr, 3));
  arr_ref<double, 2> gxsgs(cmn.gxsgs, dimension(maxstr, 3));
  arr_ref<double, 2> gysgs(cmn.gysgs, dimension(maxstr, 3));
  arr_ref<double, 2> gzsgs(cmn.gzsgs, dimension(maxstr, 3));
  arr_ref<double, 2> ftsgs(cmn.ftsgs, dimension(maxstr, 3));
  arr_ref<int, 2> k1sgs(cmn.k1sgs, dimension(maxstr, 3));
  arr_ref<int, 2> k2sgs(cmn.k2sgs, dimension(maxstr, 3));
  const int maxidl = 4001;
  arr_cref<int> itypn(cmn.itypn, dimension(maxidl));
  arr_cref<float> gxn(cmn.gxn, dimension(maxidl));
  arr_cref<float> gyn(cmn.gyn, dimension(maxidl));
  arr_cref<float> gzn(cmn.gzn, dimension(maxidl));
  arr_cref<float> ftn(cmn.ftn, dimension(maxidl));
  arr_cref<float> pxn(cmn.pxn, dimension(maxidl));
  arr_cref<float> pyn(cmn.pyn, dimension(maxidl));
  arr_cref<float> pzn(cmn.pzn, dimension(maxidl));
  arr_cref<float> een(cmn.een, dimension(maxidl));
  arr_cref<float> xmn(cmn.xmn, dimension(maxidl));
  float& bimp = cmn.bimp;
  int& iaevt = cmn.iaevt;
  int& miss = cmn.miss;
  int& ioscar = cmn.ioscar;
  float& phirp = cmn.phirp;
  arr_ref<double> xstrg0(cmn.xstrg0, dimension(maxptn));
  arr_ref<double> ystrg0(cmn.ystrg0, dimension(maxptn));
  arr_ref<int> istrg0(cmn.istrg0, dimension(maxptn));
  //
  float& aphx1 = sve.aphx1;
  float& aphx2 = sve.aphx2;
  float& b2 = sve.b2;
  float& bb = sve.bb;
  float& bbx = sve.bbx;
  float& bby = sve.bby;
  float& bmax = sve.bmax;
  float& bmin = sve.bmin;
  float& cx = sve.cx;
  float& dengy = sve.dengy;
  float& dnbp = sve.dnbp;
  float& dnbp1 = sve.dnbp1;
  float& dnbp2 = sve.dnbp2;
  float& dnbp3 = sve.dnbp3;
  float& dnbt = sve.dnbt;
  float& dnbt1 = sve.dnbt1;
  float& dnbt2 = sve.dnbt2;
  float& dnbt3 = sve.dnbt3;
  float& gs = sve.gs;
  float& gstot = sve.gstot;
  float& gstot0 = sve.gstot0;
  int& i = sve.i;
  int& i05 = sve.i05;
  int& idstr = sve.idstr;
  int& ierror = sve.ierror;
  int& ii = sve.ii;
  int& iityp = sve.iityp;
  arr_ref<int> ipcol(sve.ipcol, dimension(90000));
  int& isg = sve.isg;
  int& istr = sve.istr;
  arr_ref<int> itcol(sve.itcol, dimension(90000));
  int& itest = sve.itest;
  int& j = sve.j;
  int& jflg = sve.jflg;
  int& jjtp = sve.jjtp;
  int& jout = sve.jout;
  int& jp = sve.jp;
  int& jphard = sve.jphard;
  int& jpmini = sve.jpmini;
  int& jt = sve.jt;
  int& jthard = sve.jthard;
  int& jtmini = sve.jtmini;
  arr_ref<int> jtp(sve.jtp, dimension(3));
  int& kp = sve.kp;
  int& kp2 = sve.kp2;
  int& kt = sve.kt;
  int& kt2 = sve.kt2;
  int& lsg = sve.lsg;
  int& ncolt = sve.ncolt;
  int& nftp = sve.nftp;
  int& nhard = sve.nhard;
  int& njet = sve.njet;
  int& nlop = sve.nlop;
  int& nmini = sve.nmini;
  int& nmom = sve.nmom;
  int& npar = sve.npar;
  int& npart = sve.npart;
  int& nsbst = sve.nsbst;
  int& nsbstr = sve.nsbstr;
  int& nstrg = sve.nstrg;
  int& ntp = sve.ntp;
  float& phi = sve.phi;
  float& r = sve.r;
  float& r2 = sve.r2;
  float& rantot = sve.rantot;
  float& rnd1 = sve.rnd1;
  float& rnd2 = sve.rnd2;
  float& rnd3 = sve.rnd3;
  arr_ref<float, 2> rnip(sve.rnip, dimension(300, 300));
  float& rrb1 = sve.rrb1;
  float& rrb2 = sve.rrb2;
  arr_ref<float, 2> scip(sve.scip, dimension(300, 300));
  arr_ref<float, 2> sjip(sve.sjip, dimension(300, 300));
  float& sx = sve.sx;
  float& tt = sve.tt;
  float& ttrig = sve.ttrig;
  float& tts = sve.tts;
  float& x = sve.x;
  float& xr = sve.xr;
  float& xr1 = sve.xr1;
  float& y1 = sve.y1;
  float& y2 = sve.y2;
  float& y3 = sve.y3;
  static const char* format_210 = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))";
  static const char* format_211 = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))";
  static const char* format_243 = "(f10.3,1x,f10.3,2(1x,i5),1x,f10.3,2(1x,i5))";
  static const char* format_395 = "(3i8,f10.4,4i5)";
  //C
  //Cbz1/25/99
  //Clin-4/20/01        PARAMETER (MAXSTR = 1600)
  //Cbz1/25/99end
  //Clin-4/26/01:
  //C
  //Cbz1/31/99
  //Clin-8/2015:
  //C
  //Cbz1/31/99end
  //C
  //Cc      SAVE /HPARNT/
  //C
  //Cc      SAVE /hjcrdn/
  //Clin-7/16/03 NINT is a intrinsic fortran function, rename it to NINTHJ
  //C        COMMON/HJGLBR/NELT,NINT,NELP,NINP
  //Cc      SAVE /HJGLBR/
  //Cc      SAVE /HMAIN1/
  //Clin-4/26/01
  //C        COMMON/HMAIN2/KATT(130000,4),PATT(130000,4)
  //Cc      SAVE /HMAIN2/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /HJJET1/
  //Clin-4/2008
  //C        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
  //C     &       K2SG(900,100),PXSG(900,100),PYSG(900,100),
  //C     &       PZSG(900,100),PESG(900,100),PMSG(900,100)
  //Cc      SAVE /HJJET2/
  //Clin-4/2008:
  //C        common/xydr/rtdr(900,2)
  //Cc      SAVE /HJJET4/
  //Cc      SAVE /RNDF77/
  //C
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /LUDAT1/
  //C
  //Clin-9/29/03 changed name in order to distinguish from /prec2/
  //Ccbz11/11/98
  //C        COMMON /ARPRC/ ITYP(MAXSTR),
  //C     &     GX(MAXSTR), GY(MAXSTR), GZ(MAXSTR), FT(MAXSTR),
  //C     &     PX(MAXSTR), PY(MAXSTR), PZ(MAXSTR), EE(MAXSTR),
  //C     &     XM(MAXSTR)
  //Cc      SAVE /ARPRC/
  //Ccbz11/11/98end
  //C
  //Cbz1/25/99
  //Cc      SAVE /PARA1/
  //Cc      SAVE /prec1/
  //Cc      SAVE /prec2/
  //Cc      SAVE /ilist7/
  //Cc      SAVE /ilist8/
  //Cc      SAVE /SREC1/
  //Cc      SAVE /SREC2/
  //Cbz1/25/99end
  //C
  //Clin-2/25/00
  //Cc      SAVE /frzout/
  //Clin-4/11/01 soft:
  //Cc      SAVE /anim/
  //Clin-4/25/01 soft3:
  //Cc      SAVE /SOFT/
  //Clin-4/26/01 lepton and photon info:
  //Cc      SAVE /NOPREC/
  //Clin-6/22/01:
  //Cc      SAVE /lastt/
  //Clin-7/2011 ioscar value is needed:
  //Clin-2/2012 allow random orientation of reaction plane:
  //Clin-8/2015:
  //C
  bmax = fem::min(bmax0, hipr1(34) + hipr1(35));
  bmin = fem::min(bmin0, bmax);
  if (ihnt2(1) <= 1 && ihnt2(3) <= 1) {
    bmin = 0.0f;
    bmax = 2.5f * fem::sqrt(hipr1(31) * 0.1f / hipr1(40));
  }
  //C                        ********HIPR1(31) is in mb =0.1fm**2
  //C*******THE FOLLOWING IS TO SELECT THE COORDINATIONS OF NUCLEONS
  //C       BOTH IN PROJECTILE AND TARGET NUCLEAR( in fm)
  //C
  yp(1, 1) = 0.0f;
  yp(2, 1) = 0.0f;
  yp(3, 1) = 0.0f;
  if (ihnt2(1) <= 1) {
    goto statement_14;
  }
  FEM_DO_SAFE(kp, 1, ihnt2(1)) {
    statement_5:
    r = hirnd(cmn, 1);
    x = ranart(nseed);
    cx = 2.0f * x - 1.0f;
    sx = fem::sqrt(1.0f - cx * cx);
    //C                ********choose theta from uniform cos(theta) distr
    phi = ranart(nseed) * 2.0f * hipr1(40);
    //C                ********choose phi form uniform phi distr 0 to 2*pi
    yp(1, kp) = r * sx * fem::cos(phi);
    yp(2, kp) = r * sx * fem::sin(phi);
    yp(3, kp) = r * cx;
    if (hipr1(29) == 0.0f) {
      goto statement_10;
    }
    FEM_DO_SAFE(kp2, 1, kp - 1) {
      dnbp1 = fem::pow2((yp(1, kp) - yp(1, kp2)));
      dnbp2 = fem::pow2((yp(2, kp) - yp(2, kp2)));
      dnbp3 = fem::pow2((yp(3, kp) - yp(3, kp2)));
      dnbp = dnbp1 + dnbp2 + dnbp3;
      if (dnbp < hipr1(29) * hipr1(29)) {
        goto statement_5;
      }
      //C                        ********two neighbors cannot be closer than
      //C                                HIPR1(29)
    }
    statement_10:;
  }
  //C
  //Clin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f,
  //C     but modified [divide by 2, & x(p)=-x(n)]:
  //C     (Note: hijing1.383.f has corrected this bug in hijing1.382.f)
  if (ihnt2(1) == 2) {
    rnd1 = fem::max(ranart(nseed), 1.0e-20f);
    rnd2 = fem::max(ranart(nseed), 1.0e-20f);
    rnd3 = fem::max(ranart(nseed), 1.0e-20f);
    r = -(fem::log(rnd1) * 4.38f / 2.0f + fem::log(rnd2) * 0.85f /
      2.0f + 4.38f * 0.85f * fem::log(rnd3) / (4.38f + 0.85f));
    x = ranart(nseed);
    cx = 2.0f * x - 1.0f;
    sx = fem::sqrt(1.0f - cx * cx);
    phi = ranart(nseed) * 2.0f * hipr1(40);
    //C     R above is the relative distance between p & n in a deuteron:
    r = r / 2.f;
    yp(1, 1) = r * sx * fem::cos(phi);
    yp(2, 1) = r * sx * fem::sin(phi);
    yp(3, 1) = r * cx;
    //C     p & n has opposite coordinates in the deuteron frame:
    yp(1, 2) = -yp(1, 1);
    yp(2, 2) = -yp(2, 1);
    yp(3, 2) = -yp(3, 1);
  }
  //C
  FEM_DO_SAFE(i, 1, ihnt2(1) - 1) {
    FEM_DO_SAFE(j, i + 1, ihnt2(1)) {
      if (yp(3, i) > yp(3, j)) {
        goto statement_12;
      }
      y1 = yp(1, i);
      y2 = yp(2, i);
      y3 = yp(3, i);
      yp(1, i) = yp(1, j);
      yp(2, i) = yp(2, j);
      yp(3, i) = yp(3, j);
      yp(1, j) = y1;
      yp(2, j) = y2;
      yp(3, j) = y3;
      statement_12:;
    }
  }
  //C
  //C******************************
  statement_14:
  yt(1, 1) = 0.0f;
  yt(2, 1) = 0.0f;
  yt(3, 1) = 0.0f;
  if (ihnt2(3) <= 1) {
    goto statement_24;
  }
  FEM_DO_SAFE(kt, 1, ihnt2(3)) {
    statement_15:
    r = hirnd(cmn, 2);
    x = ranart(nseed);
    cx = 2.0f * x - 1.0f;
    sx = fem::sqrt(1.0f - cx * cx);
    //C                ********choose theta from uniform cos(theta) distr
    phi = ranart(nseed) * 2.0f * hipr1(40);
    //C                ********chose phi form uniform phi distr 0 to 2*pi
    yt(1, kt) = r * sx * fem::cos(phi);
    yt(2, kt) = r * sx * fem::sin(phi);
    yt(3, kt) = r * cx;
    if (hipr1(29) == 0.0f) {
      goto statement_20;
    }
    FEM_DO_SAFE(kt2, 1, kt - 1) {
      dnbt1 = fem::pow2((yt(1, kt) - yt(1, kt2)));
      dnbt2 = fem::pow2((yt(2, kt) - yt(2, kt2)));
      dnbt3 = fem::pow2((yt(3, kt) - yt(3, kt2)));
      dnbt = dnbt1 + dnbt2 + dnbt3;
      if (dnbt < hipr1(29) * hipr1(29)) {
        goto statement_15;
      }
      //C                        ********two neighbors cannot be closer than
      //C                                HIPR1(29)
    }
    statement_20:;
  }
  //C
  //Clin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f,
  //C     but modified [divide by 2, & x(p)=-x(n)]:
  if (ihnt2(3) == 2) {
    rnd1 = fem::max(ranart(nseed), 1.0e-20f);
    rnd2 = fem::max(ranart(nseed), 1.0e-20f);
    rnd3 = fem::max(ranart(nseed), 1.0e-20f);
    r = -(fem::log(rnd1) * 4.38f / 2.0f + fem::log(rnd2) * 0.85f /
      2.0f + 4.38f * 0.85f * fem::log(rnd3) / (4.38f + 0.85f));
    x = ranart(nseed);
    cx = 2.0f * x - 1.0f;
    sx = fem::sqrt(1.0f - cx * cx);
    phi = ranart(nseed) * 2.0f * hipr1(40);
    r = r / 2.f;
    yt(1, 1) = r * sx * fem::cos(phi);
    yt(2, 1) = r * sx * fem::sin(phi);
    yt(3, 1) = r * cx;
    yt(1, 2) = -yt(1, 1);
    yt(2, 2) = -yt(2, 1);
    yt(3, 2) = -yt(3, 1);
  }
  //C
  FEM_DO_SAFE(i, 1, ihnt2(3) - 1) {
    FEM_DO_SAFE(j, i + 1, ihnt2(3)) {
      if (yt(3, i) < yt(3, j)) {
        goto statement_22;
      }
      y1 = yt(1, i);
      y2 = yt(2, i);
      y3 = yt(3, i);
      yt(1, i) = yt(1, j);
      yt(2, i) = yt(2, j);
      yt(3, i) = yt(3, j);
      yt(1, j) = y1;
      yt(2, j) = y2;
      yt(3, j) = y3;
      statement_22:;
    }
  }
  //C
  //C********************
  statement_24:
  miss = -1;
  statement_50:
  miss++;
  //C
  //Clin-6/2009
  //C        IF(MISS.GT.50) THEN
  if (miss > cmn.maxmiss) {
    write(6, star), "infinite loop happened in  HIJING";
    FEM_STOP(0);
  }
  //C
  //Clin-4/30/01:
  itest = 0;
  //C
  natt = 0;
  jatt = 0;
  eatt = 0.0f;
  hijini(cmn);
  nlop = 0;
  //C                        ********Initialize for a new event
  statement_60:
  nt = 0;
  np = 0;
  n0 = 0;
  n01 = 0;
  n10 = 0;
  n11 = 0;
  nelt = 0;
  ninthj = 0;
  nelp = 0;
  ninp = 0;
  nsg = 0;
  ncolt = 0;
  //C
  //C****        BB IS THE ABSOLUTE VALUE OF IMPACT PARAMETER,BB**2 IS
  //C       RANDOMLY GENERATED AND ITS ORIENTATION IS RANDOMLY SET
  //C       BY THE ANGLE PHI  FOR EACH COLLISION.******************
  //C
  bb = fem::sqrt(fem::pow2(bmin) + ranart(nseed) * (fem::pow2(bmax) -
    fem::pow2(bmin)));
  //Cbz6/28/99 flow1
  //Clin-2/2012:
  phi = 0.f;
  if (cmn.iphirp == 1) {
    phi = 2.0f * hipr1(40) * ranart(nseed);
  }
  phirp = phi;
  //Cbz6/28/99 flow1 end
  bbx = bb * fem::cos(phi);
  bby = bb * fem::sin(phi);
  hint1(19) = bb;
  hint1(20) = phi;
  //C
  FEM_DO_SAFE(jp, 1, ihnt2(1)) {
    FEM_DO_SAFE(jt, 1, ihnt2(3)) {
      scip(jp, jt) = -1.0f;
      b2 = fem::pow2((yp(1, jp) + bbx - yt(1, jt))) + fem::pow2((yp(2,
        jp) + bby - yt(2, jt)));
      r2 = b2 * hipr1(40) / hipr1(31) / 0.1f;
      //C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
      rrb1 = fem::min((fem::pow2(yp(1, jp)) + fem::pow2(yp(2, jp))) /
        fem::pow2(1.2f) / fem::pow(fem::real(ihnt2(1)), 0.6666667f),
        1.0f);
      rrb2 = fem::min((fem::pow2(yt(1, jt)) + fem::pow2(yt(2, jt))) /
        fem::pow2(1.2f) / fem::pow(fem::real(ihnt2(3)), 0.6666667f),
        1.0f);
      aphx1 = hipr1(6) * 4.0f / 3.0f * (fem::pow(ihnt2(1),
        0.3333333f) - 1.0f) * fem::sqrt(1.0f - rrb1);
      aphx2 = hipr1(6) * 4.0f / 3.0f * (fem::pow(ihnt2(3),
        0.3333333f) - 1.0f) * fem::sqrt(1.0f - rrb2);
      hint1(18) = hint1(14) - aphx1 * hint1(15) - aphx2 * hint1(16) +
        aphx1 * aphx2 * hint1(17);
      if (ihpr2(14) == 0 || (ihnt2(1) == 1 && ihnt2(3) == 1)) {
        gs = 1.0f - fem::exp(-(hipr1(30) + hint1(18)) * romg(cmn,
          r2) / hipr1(31));
        rantot = ranart(nseed);
        if (rantot > gs) {
          goto statement_70;
        }
        goto statement_65;
      }
      gstot0 = 2.0f * (1.0f - fem::exp(-(hipr1(30) + hint1(18)) /
        hipr1(31) / 2.0f * romg(cmn, 0.0f)));
      r2 = r2 / gstot0;
      gs = 1.0f - fem::exp(-(hipr1(30) + hint1(18)) / hipr1(31) * romg(cmn,
        r2));
      gstot = 2.0f * (1.0f - fem::sqrt(1.0f - gs));
      rantot = ranart(nseed) * gstot0;
      if (rantot > gstot) {
        goto statement_70;
      }
      if (rantot > gs) {
        hijcsc(cmn, jp, jt);
        goto statement_70;
        //C                        ********perform elastic collisions
      }
      statement_65:
      scip(jp, jt) = r2;
      rnip(jp, jt) = rantot;
      sjip(jp, jt) = hint1(18);
      ncolt++;
      ipcol(ncolt) = jp;
      itcol(ncolt) = jt;
      statement_70:;
    }
  }
  //C                ********total number interactions proj and targ has
  //C                                suffered
  //C
  //Clin-5/22/01 write impact parameter:
  bimp = bb;
  write(6, star), "#impact parameter,nlop,ncolt=", bimp, nlop, ncolt;
  //C
  if (ncolt == 0) {
    nlop++;
    if (nlop <= 20 || (ihnt2(1) == 1 && ihnt2(3) == 1)) {
      goto statement_60;
    }
    return;
  }
  //C               ********At large impact parameter, there maybe no
  //C                       interaction at all. For NN collision
  //C                       repeat the event until interaction happens
  //C
  if (ihpr2(3) != 0) {
    nhard = 1 + fem::fint(ranart(nseed) * (ncolt - 1) + 0.5f);
    nhard = fem::min(nhard, ncolt);
    jphard = ipcol(nhard);
    jthard = itcol(nhard);
    //Clin-6/2009 ctest off:
    //C           write(99,*) IAEVT,NHARD,NCOLT,JPHARD,JTHARD
  }
  //C
  if (ihpr2(9) == 1) {
    nmini = 1 + fem::fint(ranart(nseed) * (ncolt - 1) + 0.5f);
    nmini = fem::min(nmini, ncolt);
    jpmini = ipcol(nmini);
    jtmini = itcol(nmini);
  }
  //C                ********Specifying the location of the hard and
  //C                        minijet if they are enforced by user
  //C
  FEM_DO_SAFE(jp, 1, ihnt2(1)) {
    FEM_DO_SAFE(jt, 1, ihnt2(3)) {
      if (scip(jp, jt) ==  - 1.0f) {
        goto statement_200;
      }
      nfp(jp, 11)++;
      nft(jt, 11)++;
      if (nfp(jp, 5) <= 1 && nft(jt, 5) > 1) {
        np++;
        n01++;
      }
      else if (nfp(jp, 5) > 1 && nft(jt, 5) <= 1) {
        nt++;
        n10++;
      }
      else if (nfp(jp, 5) <= 1 && nft(jt, 5) <= 1) {
        np++;
        nt++;
        n0++;
      }
      else if (nfp(jp, 5) > 1 && nft(jt, 5) > 1) {
        n11++;
      }
      jout = 0;
      nfp(jp, 10) = 0;
      nft(jt, 10) = 0;
      //C*****************************************************************
      if (ihpr2(8) == 0 && ihpr2(3) == 0) {
        goto statement_160;
      }
      //C                ********When IHPR2(8)=0 no jets are produced
      if (nfp(jp, 6) < 0 || nft(jt, 6) < 0) {
        goto statement_160;
      }
      //C                ********jets can not be produced for (JP,JT)
      //C                        because not enough energy avaible for
      //C                                JP or JT
      r2 = scip(jp, jt);
      hint1(18) = sjip(jp, jt);
      tt = romg(cmn, r2) * hint1(18) / hipr1(31);
      tts = hipr1(30) * romg(cmn, r2) / hipr1(31);
      njet = 0;
      //C
      if (ihpr2(3) != 0 && jp == jphard && jt == jthard) {
        jetini(cmn, jp, jt, 1);
        hijhrd(cmn, jp, jt, 0, jflg, 0);
        hint1(26) = hint1(47);
        hint1(27) = hint1(48);
        hint1(28) = hint1(49);
        hint1(29) = hint1(50);
        hint1(36) = hint1(67);
        hint1(37) = hint1(68);
        hint1(38) = hint1(69);
        hint1(39) = hint1(70);
        //C
        if (fem::abs(hint1(46)) > hipr1(11) && jflg == 2) {
          nfp(jp, 7) = 1;
        }
        if (fem::abs(hint1(56)) > hipr1(11) && jflg == 2) {
          nft(jt, 7) = 1;
        }
        if (fem::max(fem::abs(hint1(46)), fem::abs(hint1(
            56))) > hipr1(11) && jflg >= 3) {
          iasg(nsg, 3) = 1;
        }
        ihnt2(9) = ihnt2(14);
        ihnt2(10) = ihnt2(15);
        FEM_DO_SAFE(i05, 1, 5) {
          hint1(20 + i05) = hint1(40 + i05);
          hint1(30 + i05) = hint1(50 + i05);
        }
        //Clin-6/2009 ctest off:
        //C           write(99,*) jp,jt,IHPR2(3),HIPR1(10),njet,
        //C     1          ihnt2(9),hint1(21),hint1(22),hint1(23),
        //C     2          ihnt2(10),hint1(31),hint1(32),hint1(33)
        //C           write(99,*) ' '
        jout = 1;
        if (ihpr2(8) == 0) {
          goto statement_160;
        }
        rrb1 = fem::min((fem::pow2(yp(1, jp)) + fem::pow2(yp(2,
          jp))) / fem::pow2(1.2f) / fem::pow(fem::real(ihnt2(1)),
          0.6666667f), 1.0f);
        rrb2 = fem::min((fem::pow2(yt(1, jt)) + fem::pow2(yt(2,
          jt))) / fem::pow2(1.2f) / fem::pow(fem::real(ihnt2(3)),
          0.6666667f), 1.0f);
        aphx1 = hipr1(6) * 4.0f / 3.0f * (fem::pow(ihnt2(1),
          0.3333333f) - 1.0f) * fem::sqrt(1.0f - rrb1);
        aphx2 = hipr1(6) * 4.0f / 3.0f * (fem::pow(ihnt2(3),
          0.3333333f) - 1.0f) * fem::sqrt(1.0f - rrb2);
        hint1(65) = hint1(61) - aphx1 * hint1(62) - aphx2 * hint1(
          63) + aphx1 * aphx2 * hint1(64);
        ttrig = romg(cmn, r2) * hint1(65) / hipr1(31);
        njet = -1;
        //C                ********subtract the trigger jet from total number
        //C                        of jet production  to be done since it has
        //C                                already been produced here
        xr1 = -fem::alog(fem::exp(-ttrig) + ranart(nseed) * (1.0f -
          fem::exp(-ttrig)));
        statement_106:
        njet++;
        xr1 = xr1 - fem::alog(fem::max(ranart(nseed), 1.0e-20f));
        if (xr1 < ttrig) {
          goto statement_106;
        }
        xr = 0.0f;
        statement_107:
        njet++;
        xr = xr - fem::alog(fem::max(ranart(nseed), 1.0e-20f));
        if (xr < tt - ttrig) {
          goto statement_107;
        }
        njet = njet - 1;
        goto statement_112;
      }
      //C                ********create a hard interaction with specified P_T
      //C                                 when IHPR2(3)>0
      if (ihpr2(9) == 1 && jp == jpmini && jt == jtmini) {
        goto statement_110;
      }
      //C                ********create at least one pair of mini jets
      //C                        when IHPR2(9)=1
      //C
      //Clin-4/15/2010 changed .LT. to .LE. to avoid problem when two sides are equal;
      //C     this problem may lead to a jet production when there should be none and
      //C     crash the run; crashes at low energies were reported by P. Bhaduri.
      //C        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LT.EXP(-TT)*
      //C     &                (1.0-EXP(-TTS))) GO TO 160
      if (ihpr2(8) > 0 && rnip(jp, jt) <= fem::exp(-tt) * (1.0f -
          fem::exp(-tts))) {
        goto statement_160;
      }
      //C
      //C                ********this is the probability for no jet production
      statement_110:
      xr = -fem::alog(fem::exp(-tt) + ranart(nseed) * (1.0f - fem::exp(-tt)));
      statement_111:
      njet++;
      xr = xr - fem::alog(fem::max(ranart(nseed), 1.0e-20f));
      if (xr < tt) {
        goto statement_111;
      }
      statement_112:
      njet = fem::min(njet, ihpr2(8));
      if (ihpr2(8) < 0) {
        njet = fem::abs(ihpr2(8));
      }
      //C                ******** Determine number of mini jet production
      //C
      FEM_DO_SAFE(sve.ijet, 1, njet) {
        jetini(cmn, jp, jt, 0);
        hijhrd(cmn, jp, jt, jout, jflg, 1);
        //C                ********JFLG=1 jets valence quarks, JFLG=2 with
        //C                        gluon jet, JFLG=3 with q-qbar prod for
        //C                        (JP,JT). If JFLG=0 jets can not be produced
        //C                        this time. If JFLG=-1, error occured abandon
        //C                        this event. JOUT is the total hard scat for
        //C                        (JP,JT) up to now.
        if (jflg == 0) {
          goto statement_160;
        }
        if (jflg < 0) {
          if (ihpr2(10) != 0) {
            write(6, star), "error occured in HIJHRD";
          }
          goto statement_50;
        }
        jout++;
        if (fem::abs(hint1(46)) > hipr1(11) && jflg == 2) {
          nfp(jp, 7) = 1;
        }
        if (fem::abs(hint1(56)) > hipr1(11) && jflg == 2) {
          nft(jt, 7) = 1;
        }
        if (fem::max(fem::abs(hint1(46)), fem::abs(hint1(
            56))) > hipr1(11) && jflg >= 3) {
          iasg(nsg, 3) = 1;
        }
        //C                ******** jet with PT>HIPR1(11) will be quenched
      }
      statement_160:
      //C
      hijsft(cmn, jp, jt, jout, ierror);
      if (ierror != 0) {
        if (ihpr2(10) != 0) {
          write(6, star), "error occured in HIJSFT";
        }
        goto statement_50;
      }
      //C
      //C                ********conduct soft scattering between JP and JT
      jatt += jout;
      statement_200:;
    }
  }
  //C
  //C**************************
  //C
  //Clin-6/2009 write out initial minijet information:
  //Clin-2/2012:
  //C           call minijet_out(BB)
  minijet_out(bb, phirp);
  if (cmn.pttrig > 0 && cmn.ntrig == 0) {
    goto statement_50;
  }
  //Clin-4/2012
  //Clin-6/2009 write out initial transverse positions of initial nucleons:
  //C           write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
  FEM_DO_SAFE(jp, 1, ihnt2(1)) {
    //Clin-6/2009:
    //C           write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5)
    //Clin-2/2012:
    //C       write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5),yp(3,jp)
    //Clin-4/2012:
    //C           write(94,203) YP(1,JP)+0.5*BB*cos(phiRP),
    //C     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
    if (nfp(jp, 5) > 2) {
      ninp++;
    }
    else if (nfp(jp, 5) == 2 || nfp(jp, 5) == 1) {
      nelp++;
    }
  }
  FEM_DO_SAFE(jt, 1, ihnt2(3)) {
    //Clin-6/2009 target nucleon # has a minus sign for distinction from projectile:
    //C           write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5)
    //Clin-2/2012:
    //C       write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5),yt(3,jt)
    //Clin-4/2012:
    //C           write(94,203) YT(1,JT)-0.5*BB*cos(phiRP),
    //C     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
    if (nft(jt, 5) > 2) {
      ninthj++;
    }
    else if (nft(jt, 5) == 2 || nft(jt, 5) == 1) {
      nelt++;
    }
  }
  //C 203    format(f10.3,1x,f10.3,2(1x,I5))
  //C 203    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
  //C
  //C*******************************
  //C
  //C********perform jet quenching for jets with PT>HIPR1(11)**********
  //C
  if ((ihpr2(8) != 0 || ihpr2(3) != 0) && ihpr2(4) > 0 && ihnt2(
      1) > 1 && ihnt2(3) > 1) {
    FEM_DO_SAFE(i, 1, ihnt2(1)) {
      if (nfp(i, 7) == 1) {
        quench(cmn, i, 1);
      }
    }
    FEM_DO_SAFE(i, 1, ihnt2(3)) {
      if (nft(i, 7) == 1) {
        quench(cmn, i, 2);
      }
    }
    FEM_DO_SAFE(isg, 1, nsg) {
      if (iasg(isg, 3) == 1) {
        quench(cmn, isg, 3);
      }
    }
  }
  //C
  //Clin*****4/09/01-soft1, default way of treating strings:
  if (isoft == 1) {
    //Clin-4/16/01 allow fragmentation:
    isflag = 1;
    //C
    //Cbz1/25/99
    //C.....transfer data from HIJING to ZPC
    nsp = ihnt2(1);
    nst = ihnt2(3);
    nsi = nsg;
    istr = 0;
    npar = 0;
    FEM_DO_SAFE(i, 1, ihnt2(1)) {
      istr++;
      FEM_DO_SAFE(j, 1, npj(i)) {
        //Cbz1/27/99
        //C.....for now only consider gluon cascade
        if (kfpj(i, j) == 21) {
          //Cbz1/27/99end
          //C
          npar++;
          lstrg0(npar) = istr;
          lpart0(npar) = j;
          ityp0(npar) = kfpj(i, j);
          //Cbz6/28/99 flow1
          //Clin-7/20/01 add dble or sngl to make precisions consistent
          //C              GX0(NPAR) = YP(1, I)
          //Clin-2/2012:
          //C              GX0(NPAR) = dble(YP(1, I) + 0.5 * BB)
          gx0(npar) = fem::dble(yp(1, i) + 0.5f * bb * fem::cos(phirp));
          //Cbz6/28/99 flow1 end
          //C              GY0(NPAR) = dble(YP(2, I))
          gy0(npar) = fem::dble(yp(2, i) + 0.5f * bb * fem::sin(phirp));
          gz0(npar) = 0e0;
          ft0(npar) = 0e0;
          px0(npar) = fem::dble(pjpx(i, j));
          py0(npar) = fem::dble(pjpy(i, j));
          pz0(npar) = fem::dble(pjpz(i, j));
          xmass0(npar) = fem::dble(pjpm(i, j));
          //C              E0(NPAR) = dble(PJPE(I, J))
          e0(npar) = fem::dsqrt(fem::pow2(px0(npar)) + fem::pow2(py0(
            npar)) + fem::pow2(pz0(npar)) + fem::pow2(xmass0(npar)));
          //Clin-7/20/01-end
          //C
          //Cbz1/27/99
          //C.....end gluon selection
        }
        //Cbz1/27/99end
      }
    }
    FEM_DO_SAFE(i, 1, ihnt2(3)) {
      istr++;
      FEM_DO_SAFE(j, 1, ntj(i)) {
        //Cbz1/27/99
        //C.....for now only consider gluon cascade
        if (kftj(i, j) == 21) {
          //Cbz1/27/99end
          npar++;
          lstrg0(npar) = istr;
          lpart0(npar) = j;
          ityp0(npar) = kftj(i, j);
          //Cbz6/28/99 flow1
          //Clin-7/20/01 add dble or sngl to make precisions consistent
          //C              GX0(NPAR) = YT(1, I)
          //Clin-2/2012:
          //C              GX0(NPAR) = dble(YT(1, I) - 0.5 * BB)
          gx0(npar) = fem::dble(yt(1, i) - 0.5f * bb * fem::cos(phirp));
          //Cbz6/28/99 flow1 end
          //C              GY0(NPAR) = dble(YT(2, I))
          gy0(npar) = fem::dble(yt(2, i) - 0.5f * bb * fem::sin(phirp));
          gz0(npar) = 0e0;
          ft0(npar) = 0e0;
          px0(npar) = fem::dble(pjtx(i, j));
          py0(npar) = fem::dble(pjty(i, j));
          pz0(npar) = fem::dble(pjtz(i, j));
          xmass0(npar) = fem::dble(pjtm(i, j));
          //C              E0(NPAR) = dble(PJTE(I, J))
          e0(npar) = fem::dsqrt(fem::pow2(px0(npar)) + fem::pow2(py0(
            npar)) + fem::pow2(pz0(npar)) + fem::pow2(xmass0(npar)));
          //C
          //Cbz1/27/99
          //C.....end gluon selection
        }
        //Cbz1/27/99end
      }
    }
    FEM_DO_SAFE(i, 1, nsg) {
      istr++;
      FEM_DO_SAFE(j, 1, njsg(i)) {
        //Cbz1/27/99
        //C.....for now only consider gluon cascade
        if (k2sg(i, j) == 21) {
          //Cbz1/27/99end
          npar++;
          lstrg0(npar) = istr;
          lpart0(npar) = j;
          ityp0(npar) = k2sg(i, j);
          //Clin-7/20/01 add dble or sngl to make precisions consistent:
          gx0(npar) = 0.5e0 * fem::dble(yp(1, iasg(i, 1)) + yt(1, iasg(i, 2)));
          gy0(npar) = 0.5e0 * fem::dble(yp(2, iasg(i, 1)) + yt(2, iasg(i, 2)));
          gz0(npar) = 0e0;
          ft0(npar) = 0e0;
          px0(npar) = fem::dble(pxsg(i, j));
          py0(npar) = fem::dble(pysg(i, j));
          pz0(npar) = fem::dble(pzsg(i, j));
          xmass0(npar) = fem::dble(pmsg(i, j));
          //C              E0(NPAR) = dble(PESG(I, J))
          e0(npar) = fem::dsqrt(fem::pow2(px0(npar)) + fem::pow2(py0(
            npar)) + fem::pow2(pz0(npar)) + fem::pow2(xmass0(npar)));
          //Cbz1/27/99
          //C.....end gluon selection
        }
        //Cbz1/27/99end
      }
    }
    mul = npar;
    //C
    //Cbz2/4/99
    hjana1();
    //Cbz2/4/99end
    //C
    //Clin-6/2009:
    if (ioscar == 3) {
      write(95, star), iaevt, mul;
    }
    //C.....call ZPC for parton cascade
    zpcmn();
    //C
    //C     write out parton and wounded nucleon information to ana/zpc1.mom:
    //Clin-6/2009:
    //C        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
    write(14, format_395), iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj;
    FEM_DO_SAFE(i, 1, mul) {
      //Cc           WRITE (14, 411) PX5(I), PY5(I), PZ5(I), ITYP5(I),
      //C     &        XMASS5(I), E5(I)
      if (fem::dmax1(fem::abs(gx5(i)), fem::abs(gy5(i)), fem::abs(gz5(i)),
          fem::abs(ft5(i))) < 9999) {
        write(14, format_210), ityp5(i), px5(i), py5(i), pz5(i),
          xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i);
      }
      else {
        //C     change format for large numbers:
        write(14, format_211), ityp5(i), px5(i), py5(i), pz5(i),
          xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i);
      }
      //C
    }
    //C
    //Clin-4/09/01:
    itest++;
    //C 411    FORMAT(1X, 3F10.3, I6, 2F10.3)
    //Cbz3/19/99 end
    //C
    //Clin-5/2009 ctest off:
    //C        call frztm(1,1)
    //C
    //C.....transfer data back from ZPC to HIJING
    FEM_DO_SAFE(i, 1, mul) {
      if (lstrg1(i) <= nsp) {
        nstrg = lstrg1(i);
        npart = lpart1(i);
        kfpj(nstrg, npart) = ityp5(i);
        //Clin-7/20/01 add dble or sngl to make precisions consistent
        pjpx(nstrg, npart) = fem::sngl(px5(i));
        pjpy(nstrg, npart) = fem::sngl(py5(i));
        pjpz(nstrg, npart) = fem::sngl(pz5(i));
        pjpe(nstrg, npart) = fem::sngl(e5(i));
        pjpm(nstrg, npart) = fem::sngl(xmass5(i));
      }
      else if (lstrg1(i) <= nsp + nst) {
        nstrg = lstrg1(i) - nsp;
        npart = lpart1(i);
        kftj(nstrg, npart) = ityp5(i);
        pjtx(nstrg, npart) = fem::sngl(px5(i));
        pjty(nstrg, npart) = fem::sngl(py5(i));
        pjtz(nstrg, npart) = fem::sngl(pz5(i));
        pjte(nstrg, npart) = fem::sngl(e5(i));
        pjtm(nstrg, npart) = fem::sngl(xmass5(i));
      }
      else {
        nstrg = lstrg1(i) - nsp - nst;
        npart = lpart1(i);
        k2sg(nstrg, npart) = ityp5(i);
        pxsg(nstrg, npart) = fem::sngl(px5(i));
        pysg(nstrg, npart) = fem::sngl(py5(i));
        pzsg(nstrg, npart) = fem::sngl(pz5(i));
        pesg(nstrg, npart) = fem::sngl(e5(i));
        pmsg(nstrg, npart) = fem::sngl(xmass5(i));
      }
    }
    //Cbz1/25/99end
    //C
    //Cbz2/4/99
    hjana2();
    //Cbz2/4/99end
    //C
    //Clin*****4/09/01-soft2, put q+dq+X in strings into ZPC:
  }
  else if (isoft == 2) {
    nsp = ihnt2(1);
    nst = ihnt2(3);
    //Clin-4/27/01:
    nsi = nsg;
    npar = 0;
    istr = 0;
    //C
    //Clin  No fragmentation to hadrons, only on parton level,
    //C     and transfer minijet and string data from HIJING to ZPC:
    mstj(1) = 0;
    //Clin-4/12/01 forbid soft radiation before ZPC to avoid small-mass strings,
    //C     and forbid jet order reversal before ZPC to avoid unphysical flavors:
    ihpr2(1) = 0;
    isflag = 0;
    //C
    if (ihpr2(20) != 0) {
      FEM_DO_SAFE(ntp, 1, 2) {
        FEM_DO_SAFE(jjtp, 1, ihnt2(2 * ntp - 1)) {
          istr++;
          //C change: do gluon kink only once: either here or in fragmentation.
          hijfrg(cmn, jjtp, ntp, ierror);
          //C                 call lulist(1)
          if (ntp == 1) {
            //C 354                continue
            npj(jjtp) = fem::max0(n - 2, 0);
            //C
            //Clin-4/12/01:                    NPJ(jjtp)=MAX0(ipartn-2,0)
          }
          else {
            //C 355                continue
            ntj(jjtp) = fem::max0(n - 2, 0);
            //Clin-4/12/01:                    NTJ(jjtp)=MAX0(ipartn-2,0)
          }
          //C
          FEM_DO_SAFE(ii, 1, n) {
            npar++;
            lstrg0(npar) = istr;
            lpart0(npar) = ii;
            ityp0(npar) = k(ii, 2);
            gz0(npar) = 0e0;
            ft0(npar) = 0e0;
            //Clin-7/20/01 add dble or sngl to make precisions consistent
            px0(npar) = fem::dble(p(ii, 1));
            py0(npar) = fem::dble(p(ii, 2));
            pz0(npar) = fem::dble(p(ii, 3));
            xmass0(npar) = fem::dble(p(ii, 5));
            //C                 E0(NPAR) = dble(P(II,4))
            e0(npar) = fem::dsqrt(fem::pow2(px0(npar)) + fem::pow2(
              py0(npar)) + fem::pow2(pz0(npar)) + fem::pow2(xmass0(
              npar)));
            if (ntp == 1) {
              //Clin-7/20/01 add dble or sngl to make precisions consistent
              //Clin-2/2012:
              //C                    GX0(NPAR) = dble(YP(1, jjtp)+0.5 * BB)
              //C                    GY0(NPAR) = dble(YP(2, jjtp))
              gx0(npar) = fem::dble(yp(1, jjtp) + 0.5f * bb * fem::cos(phirp));
              gy0(npar) = fem::dble(yp(2, jjtp) + 0.5f * bb * fem::sin(phirp));
              //C
              iityp = ityp0(npar);
              nstrg = lstrg0(npar);
              if (iityp == 2112 || iityp == 2212) {
              }
              else if ((iityp == 1 || iityp == 2) && (ii == 1 || ii == n)) {
                pp(nstrg, 6) = fem::sngl(px0(npar));
                pp(nstrg, 7) = fem::sngl(py0(npar));
                pp(nstrg, 14) = fem::sngl(xmass0(npar));
              }
              else if ((iityp == 1103 || iityp == 2101 ||
                iityp == 2103 || iityp == 2203.f || iityp == 3101 ||
                iityp == 3103.f || iityp == 3201 || iityp == 3203 ||
                iityp == 3303) && (ii == 1 || ii == n)) {
                pp(nstrg, 8) = fem::sngl(px0(npar));
                pp(nstrg, 9) = fem::sngl(py0(npar));
                pp(nstrg, 15) = fem::sngl(xmass0(npar));
              }
              else {
                npart = lpart0(npar) - 1;
                kfpj(nstrg, npart) = ityp0(npar);
                pjpx(nstrg, npart) = fem::sngl(px0(npar));
                pjpy(nstrg, npart) = fem::sngl(py0(npar));
                pjpz(nstrg, npart) = fem::sngl(pz0(npar));
                pjpe(nstrg, npart) = fem::sngl(e0(npar));
                pjpm(nstrg, npart) = fem::sngl(xmass0(npar));
              }
            }
            else {
              //Clin-2/2012:
              //C                    GX0(NPAR) = dble(YT(1, jjtp)-0.5 * BB)
              //C                    GY0(NPAR) = dble(YT(2, jjtp))
              gx0(npar) = fem::dble(yt(1, jjtp) - 0.5f * bb * fem::cos(phirp));
              gy0(npar) = fem::dble(yt(2, jjtp) - 0.5f * bb * fem::sin(phirp));
              iityp = ityp0(npar);
              nstrg = lstrg0(npar) - nsp;
              if (iityp == 2112 || iityp == 2212) {
              }
              else if ((iityp == 1 || iityp == 2) && (ii == 1 || ii == n)) {
                pt(nstrg, 6) = fem::sngl(px0(npar));
                pt(nstrg, 7) = fem::sngl(py0(npar));
                pt(nstrg, 14) = fem::sngl(xmass0(npar));
              }
              else if ((iityp == 1103 || iityp == 2101 ||
                iityp == 2103 || iityp == 2203.f || iityp == 3101 ||
                iityp == 3103.f || iityp == 3201 || iityp == 3203 ||
                iityp == 3303) && (ii == 1 || ii == n)) {
                pt(nstrg, 8) = fem::sngl(px0(npar));
                pt(nstrg, 9) = fem::sngl(py0(npar));
                pt(nstrg, 15) = fem::sngl(xmass0(npar));
              }
              else {
                npart = lpart0(npar) - 1;
                kftj(nstrg, npart) = ityp0(npar);
                pjtx(nstrg, npart) = fem::sngl(px0(npar));
                pjty(nstrg, npart) = fem::sngl(py0(npar));
                pjtz(nstrg, npart) = fem::sngl(pz0(npar));
                pjte(nstrg, npart) = fem::sngl(e0(npar));
                pjtm(nstrg, npart) = fem::sngl(xmass0(npar));
              }
            }
          }
        }
      }
      FEM_DO_SAFE(isg, 1, nsg) {
        istr++;
        hijfrg(cmn, isg, 3, ierror);
        //C              call lulist(2)
        //C
        njsg(isg) = n;
        //C
        FEM_DO_SAFE(ii, 1, n) {
          npar++;
          lstrg0(npar) = istr;
          lpart0(npar) = ii;
          ityp0(npar) = k(ii, 2);
          gx0(npar) = 0.5e0 * fem::dble(yp(1, iasg(isg, 1)) + yt(1,
            iasg(isg, 2)));
          gy0(npar) = 0.5e0 * fem::dble(yp(2, iasg(isg, 1)) + yt(2,
            iasg(isg, 2)));
          gz0(npar) = 0e0;
          ft0(npar) = 0e0;
          px0(npar) = fem::dble(p(ii, 1));
          py0(npar) = fem::dble(p(ii, 2));
          pz0(npar) = fem::dble(p(ii, 3));
          xmass0(npar) = fem::dble(p(ii, 5));
          //C                 E0(NPAR) = dble(P(II,4))
          e0(npar) = fem::dsqrt(fem::pow2(px0(npar)) + fem::pow2(py0(
            npar)) + fem::pow2(pz0(npar)) + fem::pow2(xmass0(npar)));
        }
      }
    }
    //C
    mul = npar;
    //Cbz2/4/99
    hjana1();
    //Cbz2/4/99end
    //Clin-6/2009:
    if (ioscar == 3) {
      write(95, star), iaevt, mul;
    }
    //C.....call ZPC for parton cascade
    zpcmn();
    //Cbz3/19/99
    //Clin-6/2009:
    //C        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
    write(14, format_395), iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj;
    itest++;
    //C
    FEM_DO_SAFE(i, 1, mul) {
      //C           WRITE (14, 311) PX5(I), PY5(I), PZ5(I), ITYP5(I),
      //C     &        XMASS5(I), E5(I)
      //Clin-4/2012 write parton freeze-out position in zpc.dat for this test scenario:
      //C           WRITE (14, 312) PX5(I), PY5(I), PZ5(I), ITYP5(I),
      //C     &        XMASS5(I), E5(I),LSTRG1(I), LPART1(I)
      if (fem::dmax1(fem::abs(gx5(i)), fem::abs(gy5(i)), fem::abs(gz5(i)),
          fem::abs(ft5(i))) < 9999) {
        write(14, format_210), ityp5(i), px5(i), py5(i), pz5(i),
          xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i);
      }
      else {
        write(14, format_211), ityp5(i), px5(i), py5(i), pz5(i),
          xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i);
      }
      //C
    }
    //C 311    FORMAT(1X, 3F10.4, I6, 2F10.4)
    //C 312    FORMAT(1X, 3F10.3, I6, 2F10.3,1X,I6,1X,I3)
    //Cbz3/19/99 end
    //C
    //Clin-5/2009 ctest off:
    //C        call frztm(1,1)
    //C
    //Clin-4/13/01 initialize four momenta and invariant mass of strings after ZPC:
    FEM_DO_SAFE(nmom, 1, 5) {
      FEM_DO_SAFE(nstrg, 1, nsp) {
        pp(nstrg, nmom) = 0.f;
      }
      FEM_DO_SAFE(nstrg, 1, nst) {
        pt(nstrg, nmom) = 0.f;
      }
    }
    //Clin-4/13/01-end
    //C
    FEM_DO_SAFE(i, 1, mul) {
      iityp = ityp5(i);
      if (lstrg1(i) <= nsp) {
        nstrg = lstrg1(i);
        //C     nucleons without interactions:
        if (iityp == 2112 || iityp == 2212) {
          //Clin-7/20/01 add dble or sngl to make precisions consistent
          pp(nstrg, 1) = fem::sngl(px5(i));
          pp(nstrg, 2) = fem::sngl(py5(i));
          pp(nstrg, 3) = fem::sngl(pz5(i));
          pp(nstrg, 4) = fem::sngl(e5(i));
          pp(nstrg, 5) = fem::sngl(xmass5(i));
          //C     valence quark:
        }
        else if ((iityp == 1 || iityp == 2) && (lpart1(i) == 1 ||
          lpart1(i) == (npj(nstrg) + 2))) {
          pp(nstrg, 6) = fem::sngl(px5(i));
          pp(nstrg, 7) = fem::sngl(py5(i));
          pp(nstrg, 14) = fem::sngl(xmass5(i));
          pp(nstrg, 1) += fem::sngl(px5(i));
          pp(nstrg, 2) += fem::sngl(py5(i));
          pp(nstrg, 3) += fem::sngl(pz5(i));
          pp(nstrg, 4) += fem::sngl(e5(i));
          pp(nstrg, 5) = fem::sqrt(fem::pow2(pp(nstrg, 4)) -
            fem::pow2(pp(nstrg, 1)) - fem::pow2(pp(nstrg, 2)) -
            fem::pow2(pp(nstrg, 3)));
          //C     diquark:
        }
        else if ((iityp == 1103 || iityp == 2101 || iityp == 2103 ||
          iityp == 2203.f || iityp == 3101 || iityp == 3103.f ||
          iityp == 3201 || iityp == 3203 || iityp == 3303) && (lpart1(
          i) == 1 || lpart1(i) == (npj(nstrg) + 2))) {
          pp(nstrg, 8) = fem::sngl(px5(i));
          pp(nstrg, 9) = fem::sngl(py5(i));
          pp(nstrg, 15) = fem::sngl(xmass5(i));
          pp(nstrg, 1) += fem::sngl(px5(i));
          pp(nstrg, 2) += fem::sngl(py5(i));
          pp(nstrg, 3) += fem::sngl(pz5(i));
          pp(nstrg, 4) += fem::sngl(e5(i));
          pp(nstrg, 5) = fem::sqrt(fem::pow2(pp(nstrg, 4)) -
            fem::pow2(pp(nstrg, 1)) - fem::pow2(pp(nstrg, 2)) -
            fem::pow2(pp(nstrg, 3)));
          //C     partons in projectile or target strings:
        }
        else {
          npart = lpart1(i) - 1;
          kfpj(nstrg, npart) = ityp5(i);
          pjpx(nstrg, npart) = fem::sngl(px5(i));
          pjpy(nstrg, npart) = fem::sngl(py5(i));
          pjpz(nstrg, npart) = fem::sngl(pz5(i));
          pjpe(nstrg, npart) = fem::sngl(e5(i));
          pjpm(nstrg, npart) = fem::sngl(xmass5(i));
        }
      }
      else if (lstrg1(i) <= nsp + nst) {
        nstrg = lstrg1(i) - nsp;
        if (iityp == 2112 || iityp == 2212) {
          pt(nstrg, 1) = fem::sngl(px5(i));
          pt(nstrg, 2) = fem::sngl(py5(i));
          pt(nstrg, 3) = fem::sngl(pz5(i));
          pt(nstrg, 4) = fem::sngl(e5(i));
          pt(nstrg, 5) = fem::sngl(xmass5(i));
        }
        else if ((iityp == 1 || iityp == 2) && (lpart1(i) == 1 ||
          lpart1(i) == (ntj(nstrg) + 2))) {
          pt(nstrg, 6) = fem::sngl(px5(i));
          pt(nstrg, 7) = fem::sngl(py5(i));
          pt(nstrg, 14) = fem::sngl(xmass5(i));
          pt(nstrg, 1) += fem::sngl(px5(i));
          pt(nstrg, 2) += fem::sngl(py5(i));
          pt(nstrg, 3) += fem::sngl(pz5(i));
          pt(nstrg, 4) += fem::sngl(e5(i));
          pt(nstrg, 5) = fem::sqrt(fem::pow2(pt(nstrg, 4)) -
            fem::pow2(pt(nstrg, 1)) - fem::pow2(pt(nstrg, 2)) -
            fem::pow2(pt(nstrg, 3)));
        }
        else if ((iityp == 1103 || iityp == 2101 || iityp == 2103 ||
          iityp == 2203.f || iityp == 3101 || iityp == 3103.f ||
          iityp == 3201 || iityp == 3203 || iityp == 3303) && (lpart1(
          i) == 1 || lpart1(i) == (ntj(nstrg) + 2))) {
          pt(nstrg, 8) = fem::sngl(px5(i));
          pt(nstrg, 9) = fem::sngl(py5(i));
          pt(nstrg, 15) = fem::sngl(xmass5(i));
          pt(nstrg, 1) += fem::sngl(px5(i));
          pt(nstrg, 2) += fem::sngl(py5(i));
          pt(nstrg, 3) += fem::sngl(pz5(i));
          pt(nstrg, 4) += fem::sngl(e5(i));
          pt(nstrg, 5) = fem::sqrt(fem::pow2(pt(nstrg, 4)) -
            fem::pow2(pt(nstrg, 1)) - fem::pow2(pt(nstrg, 2)) -
            fem::pow2(pt(nstrg, 3)));
        }
        else {
          npart = lpart1(i) - 1;
          kftj(nstrg, npart) = ityp5(i);
          pjtx(nstrg, npart) = fem::sngl(px5(i));
          pjty(nstrg, npart) = fem::sngl(py5(i));
          pjtz(nstrg, npart) = fem::sngl(pz5(i));
          pjte(nstrg, npart) = fem::sngl(e5(i));
          pjtm(nstrg, npart) = fem::sngl(xmass5(i));
        }
      }
      else {
        nstrg = lstrg1(i) - nsp - nst;
        npart = lpart1(i);
        k2sg(nstrg, npart) = ityp5(i);
        pxsg(nstrg, npart) = fem::sngl(px5(i));
        pysg(nstrg, npart) = fem::sngl(py5(i));
        pzsg(nstrg, npart) = fem::sngl(pz5(i));
        pesg(nstrg, npart) = fem::sngl(e5(i));
        pmsg(nstrg, npart) = fem::sngl(xmass5(i));
      }
    }
    //Cbz1/25/99end
    //C
    //Clin-4/09/01  turn on fragmentation with soft radiation
    //C     and jet order reversal to form hadrons after ZPC:
    mstj(1) = 1;
    ihpr2(1) = 1;
    isflag = 1;
    //Clin-4/13/01 allow small mass strings (D=1.5GeV):
    hipr1(1) = 0.94f;
    //C
    //Cbz2/4/99
    hjana2();
    //Cbz2/4/99end
    //C
    //Clin-4/19/01-soft3, fragment strings, then convert hadrons to partons
    //C     and input to ZPC:
  }
  else if (isoft == 3 || isoft == 4 || isoft == 5) {
    //Clin-4/24/01 normal fragmentation first:
    isflag = 0;
    //C        write(99,*) 'IAEVT,NSG,NDR=',IAEVT,NSG,NDR
    //C
    if (ihpr2(20) != 0) {
      FEM_DO_SAFE(isg, 1, nsg) {
        hijfrg(cmn, isg, 3, ierror);
        //C
        nsbst = 1;
        idstr = 92;
        if (ihpr2(21) == 0) {
          luedit(2);
        }
        else {
          statement_551:
          nsbst++;
          if (k(nsbst, 2) < 91 || k(nsbst, 2) > 93) {
            goto statement_551;
          }
          idstr = k(nsbst, 2);
          nsbst++;
        }
        //C
        if (frame == "LAB") {
          hboost(cmn);
        }
        //C                ******** boost back to lab frame(if it was in)
        //C
        nsbstr = 0;
        FEM_DO_SAFE(i, nsbst, n) {
          if (k(i, 2) == idstr) {
            nsbstr++;
            goto statement_560;
          }
          k(i, 4) = nsbstr;
          natt++;
          katt(natt, 1) = k(i, 2);
          katt(natt, 2) = 20;
          katt(natt, 4) = k(i, 1);
          //C     from Yasushi, to avoid violation of array limits:
          //C                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
          //Clin-4/2008 to avoid out-of-bound error in K():
          //C                   IF(K(I,3).EQ.0 .OR.
          //C     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
          //C                      KATT(NATT,3)=0
          if (k(i, 3) == 0) {
            katt(natt, 3) = 0;
          }
          else if (k(i, 3) != 0 && k(k(i, 3), 2) == idstr) {
            katt(natt, 3) = 0;
            //Clin-4/2008-end
          }
          else {
            katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i, 3), 4);
          }
          //C
          //C       ****** identify the mother particle
          patt(natt, 1) = p(i, 1);
          patt(natt, 2) = p(i, 2);
          patt(natt, 3) = p(i, 3);
          patt(natt, 4) = p(i, 4);
          eatt += p(i, 4);
          gxar(natt) = 0.5f * (yp(1, iasg(isg, 1)) + yt(1, iasg(isg, 2)));
          gyar(natt) = 0.5f * (yp(2, iasg(isg, 1)) + yt(2, iasg(isg, 2)));
          gzar(natt) = 0.f;
          ftar(natt) = 0.f;
          itypar(natt) = k(i, 2);
          pxar(natt) = p(i, 1);
          pyar(natt) = p(i, 2);
          pzar(natt) = p(i, 3);
          pear(natt) = p(i, 4);
          xmar(natt) = p(i, 5);
          //Clin-8/2015: record hadron information, to be used for its constituent partons:
          xstrg0(natt) = fem::dble(gxar(natt));
          ystrg0(natt) = fem::dble(gyar(natt));
          istrg0(natt) = isg;
          //C                   write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT),
          //C     1                  K(I,2),P(I, 1),P(I, 2),P(I, 3)
          //Cbz11/11/98end
          //C
          statement_560:;
        }
      }
      //C                ********Fragment the q-qbar jets systems *****
      //C
      jtp(1) = ihnt2(1);
      jtp(2) = ihnt2(3);
      FEM_DO_SAFE(ntp, 1, 2) {
        FEM_DO_SAFE(jjtp, 1, jtp(ntp)) {
          hijfrg(cmn, jjtp, ntp, ierror);
          //C
          nsbst = 1;
          idstr = 92;
          if (ihpr2(21) == 0) {
            luedit(2);
          }
          else {
            statement_581:
            nsbst++;
            if (k(nsbst, 2) < 91 || k(nsbst, 2) > 93) {
              goto statement_581;
            }
            idstr = k(nsbst, 2);
            nsbst++;
          }
          if (frame == "LAB") {
            hboost(cmn);
          }
          //C                ******** boost back to lab frame(if it was in)
          //C
          nftp = nfp(jjtp, 5);
          if (ntp == 2) {
            nftp = 10 + nft(jjtp, 5);
          }
          nsbstr = 0;
          FEM_DO_SAFE(i, nsbst, n) {
            if (k(i, 2) == idstr) {
              nsbstr++;
              goto statement_590;
            }
            k(i, 4) = nsbstr;
            natt++;
            katt(natt, 1) = k(i, 2);
            katt(natt, 2) = nftp;
            katt(natt, 4) = k(i, 1);
            //C                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
            //Clin-4/2008
            //C                   IF(K(I,3).EQ.0 .OR.
            //C     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
            //C                      KATT(NATT,3)=0
            if (k(i, 3) == 0) {
              katt(natt, 3) = 0;
            }
            else if (k(i, 3) != 0 && k(k(i, 3), 2) == idstr) {
              katt(natt, 3) = 0;
              //Clin-4/2008-end
            }
            else {
              katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i, 3), 4);
            }
            //C
            //C       ****** identify the mother particle
            patt(natt, 1) = p(i, 1);
            patt(natt, 2) = p(i, 2);
            patt(natt, 3) = p(i, 3);
            patt(natt, 4) = p(i, 4);
            eatt += p(i, 4);
            if (ntp == 1) {
              //Clin-2/2012:
              //C                      GXAR(NATT) = YP(1, jjtp)+0.5 * BB
              //C                      GYAR(NATT) = YP(2, jjtp)
              gxar(natt) = yp(1, jjtp) + 0.5f * bb * fem::cos(phirp);
              gyar(natt) = yp(2, jjtp) + 0.5f * bb * fem::sin(phirp);
              //C
            }
            else {
              //Clin-2/2012:
              //C                      GXAR(NATT) = YT(1, jjtp)-0.5 * BB
              //C                      GYAR(NATT) = YT(2, jjtp)
              gxar(natt) = yt(1, jjtp) - 0.5f * bb * fem::cos(phirp);
              gyar(natt) = yt(2, jjtp) - 0.5f * bb * fem::sin(phirp);
            }
            gzar(natt) = 0.f;
            ftar(natt) = 0.f;
            itypar(natt) = k(i, 2);
            pxar(natt) = p(i, 1);
            pyar(natt) = p(i, 2);
            pzar(natt) = p(i, 3);
            pear(natt) = p(i, 4);
            xmar(natt) = p(i, 5);
            //Clin-8/2015: record hadron information, to be used for its constituent partons:
            xstrg0(natt) = fem::dble(gxar(natt));
            ystrg0(natt) = fem::dble(gyar(natt));
            //C     String ID is separated for projectile/target strings:
            istrg0(natt) = ntp * 10000 + jjtp;
            //C              if(N.eq.nsbst.and.(K(I,2).eq.2112.or.K(I,2).eq.2212)) then
            //C                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
            //C     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3),'spectator'
            //C                   else
            //C                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
            //C     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3)
            //C                   endif
            //Cbz11/11/98end
            //C
            statement_590:;
          }
        }
      }
      //C     ********Fragment the q-qq related string systems
    }
    //Clin-4/2008 check for zero NDR value:
    if (ndr >= 1) {
      //C
      FEM_DO_SAFE(i, 1, ndr) {
        natt++;
        katt(natt, 1) = kfdr(i);
        katt(natt, 2) = 40;
        katt(natt, 3) = 0;
        patt(natt, 1) = pdr(i, 1);
        patt(natt, 2) = pdr(i, 2);
        patt(natt, 3) = pdr(i, 3);
        patt(natt, 4) = pdr(i, 4);
        eatt += pdr(i, 4);
        //Clin-11/11/03     set direct photons positions and time at formation:
        gxar(natt) = rtdr(i, 1);
        gyar(natt) = rtdr(i, 2);
        gzar(natt) = 0.f;
        ftar(natt) = 0.f;
        itypar(natt) = katt(natt, 1);
        pxar(natt) = patt(natt, 1);
        pyar(natt) = patt(natt, 2);
        pzar(natt) = patt(natt, 3);
        pear(natt) = patt(natt, 4);
        xmar(natt) = pdr(i, 5);
      }
      //Clin-4/2008:
    }
    //Clin-6/2009
    embedhighpt();
    //C
    hjana1();
    //C
    //Clin-4/19/01 convert hadrons to partons for ZPC (with GX0 given):
    htop();
    //C
    //Clin-7/03/01 move up, used in zpstrg (otherwise not set and incorrect):
    nsp = 0;
    nst = 0;
    nsg = natt;
    nsi = nsg;
    //Clin-7/03/01-end
    //C
    //Clin-6/2009:
    if (ioscar == 3) {
      write(95, star), iaevt, mul;
    }
    //C
    //C.....call ZPC for parton cascade
    zpcmn();
    //Clin-6/2009:
    //C        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
    write(14, format_395), iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj;
    itest++;
    //C
    FEM_DO_SAFE(i, 1, mul) {
      //C           WRITE (14, 511) PX5(I), PY5(I), PZ5(I), ITYP5(I),
      //C     &        XMASS5(I), E5(I)
      //Clin-4/2012 write parton freeze-out position in zpc.dat
      //C     for string melting version:
      //C           WRITE (14, 512) ITYP5(I), PX5(I), PY5(I), PZ5(I),
      //C     &        XMASS5(I), LSTRG1(I), LPART1(I), FT5(I)
      if (fem::dmax1(fem::abs(gx5(i)), fem::abs(gy5(i)), fem::abs(gz5(i)),
          fem::abs(ft5(i))) < 9999) {
        write(14, format_210), ityp5(i), px5(i), py5(i), pz5(i),
          xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i);
      }
      else {
        write(14, format_211), ityp5(i), px5(i), py5(i), pz5(i),
          xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i);
      }
      //C
    }
    //C 511    FORMAT(1X, 3F10.4, I6, 2F10.4)
    //C 512    FORMAT(I6,4(1X,F10.3),1X,I6,1X,I3,1X,F10.3)
    //C 513    FORMAT(1X, 4F10.4)
    //C
    //Clin-5/2009 ctest off:
    //C        call frztm(1,1)
    //C
    //Clin  save data after ZPC for fragmentation purpose:
    //C.....transfer data back from ZPC to HIJING
    FEM_DO_SAFE(i, 1, maxstr) {
      FEM_DO_SAFE(j, 1, 3) {
        k1sgs(i, j) = 0;
        k2sgs(i, j) = 0;
        pxsgs(i, j) = 0e0;
        pysgs(i, j) = 0e0;
        pzsgs(i, j) = 0e0;
        pesgs(i, j) = 0e0;
        pmsgs(i, j) = 0e0;
        gxsgs(i, j) = 0e0;
        gysgs(i, j) = 0e0;
        gzsgs(i, j) = 0e0;
        ftsgs(i, j) = 0e0;
      }
    }
    FEM_DO_SAFE(i, 1, mul) {
      iityp = ityp5(i);
      nstrg = lstrg1(i);
      npart = lpart1(i);
      k2sgs(nstrg, npart) = ityp5(i);
      pxsgs(nstrg, npart) = px5(i);
      pysgs(nstrg, npart) = py5(i);
      pzsgs(nstrg, npart) = pz5(i);
      pmsgs(nstrg, npart) = xmass5(i);
      //Clin-7/20/01 E5(I) does no include the finite parton mass XMASS5(I),
      //C     so define it anew:
      //C           PESGS(NSTRG, NPART) = E5(I)
      //C           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0)
      //C     1          write(91,*) 'a',PX5(i),PY5(i),XMASS5(i),PZ5(i),E5(i)
      e5(i) = fem::dsqrt(fem::pow2(px5(i)) + fem::pow2(py5(i)) +
        fem::pow2(pz5(i)) + fem::pow2(xmass5(i)));
      pesgs(nstrg, npart) = e5(i);
      //C           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0)
      //C     1          write(91,*) 'b: new E5(I)=',E5(i)
      //Clin-7/20/01-end
      gxsgs(nstrg, npart) = gx5(i);
      gysgs(nstrg, npart) = gy5(i);
      gzsgs(nstrg, npart) = gz5(i);
      ftsgs(nstrg, npart) = ft5(i);
    }
    hjana2();
    //C
    //Clin-4/19/01-end
    //C
  }
  //Clin-4/09/01-end
  //C
  //C**************fragment all the string systems in the following*****
  //C
  //C********nsbst is where particle information starts
  //C********nsbstR+1 is the number of strings in fragmentation
  //C********the number of strings before a line is stored in K(I,4)
  //C********IDSTR is id number of the string system (91,92 or 93)
  //C
  //Clin-4/30/01 convert partons to hadrons after ZPC:
  if (isoft == 3 || isoft == 4 || isoft == 5) {
    natt = 0;
    eatt = 0.f;
    ptoh();
    FEM_DO_SAFE(i, 1, cmn.nnozpc) {
      natt++;
      katt(natt, 1) = itypn(i);
      patt(natt, 1) = pxn(i);
      patt(natt, 2) = pyn(i);
      patt(natt, 3) = pzn(i);
      patt(natt, 4) = een(i);
      eatt += een(i);
      gxar(natt) = gxn(i);
      gyar(natt) = gyn(i);
      gzar(natt) = gzn(i);
      ftar(natt) = ftn(i);
      itypar(natt) = itypn(i);
      pxar(natt) = pxn(i);
      pyar(natt) = pyn(i);
      pzar(natt) = pzn(i);
      pear(natt) = een(i);
      xmar(natt) = xmn(i);
    }
    goto statement_565;
  }
  //Clin-4/30/01-end
  if (ihpr2(20) != 0) {
    FEM_DO_SAFE(isg, 1, nsg) {
      hijfrg(cmn, isg, 3, ierror);
      if (mstu(24) != 0 || ierror > 0) {
        mstu(24) = 0;
        mstu(28) = 0;
        if (ihpr2(10) != 0) {
          //C                      call lulist(2)
          write(6, star), "error occured ISG, repeat the event";
          write(6, star), isg;
          //C
        }
        goto statement_50;
      }
      //C                        ********Check errors
      //C
      nsbst = 1;
      idstr = 92;
      if (ihpr2(21) == 0) {
        luedit(2);
      }
      else {
        statement_351:
        nsbst++;
        if (k(nsbst, 2) < 91 || k(nsbst, 2) > 93) {
          goto statement_351;
        }
        idstr = k(nsbst, 2);
        nsbst++;
      }
      //C
      if (frame == "LAB") {
        hboost(cmn);
      }
      //C                ******** boost back to lab frame(if it was in)
      //C
      nsbstr = 0;
      FEM_DO_SAFE(i, nsbst, n) {
        if (k(i, 2) == idstr) {
          nsbstr++;
          goto statement_360;
        }
        k(i, 4) = nsbstr;
        natt++;
        katt(natt, 1) = k(i, 2);
        katt(natt, 2) = 20;
        katt(natt, 4) = k(i, 1);
        //C                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
        //Clin-4/2008:
        //C                   IF(K(I,3).EQ.0 .OR.
        //C     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
        //C                      KATT(NATT,3)=0
        if (k(i, 3) == 0) {
          katt(natt, 3) = 0;
        }
        else if (k(i, 3) != 0 && k(k(i, 3), 2) == idstr) {
          katt(natt, 3) = 0;
          //Clin-4/2008-end
        }
        else {
          katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i, 3), 4);
        }
        //C
        //C       ****** identify the mother particle
        patt(natt, 1) = p(i, 1);
        patt(natt, 2) = p(i, 2);
        patt(natt, 3) = p(i, 3);
        patt(natt, 4) = p(i, 4);
        eatt += p(i, 4);
        //C
        //Cbz11/11/98
        //Cbz1/25/99
        //C                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
        //C     &                YT(1, IASG(ISG, 2)))
        //C                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
        //C     &                YT(2, IASG(ISG, 2)))
        lsg = nsp + nst + isg;
        gxar(natt) = fem::sngl(zt1(lsg));
        gyar(natt) = fem::sngl(zt2(lsg));
        gzar(natt) = fem::sngl(zt3(lsg));
        ftar(natt) = fem::sngl(ataui(lsg));
        //Cbz1/25/99end
        itypar(natt) = k(i, 2);
        pxar(natt) = p(i, 1);
        pyar(natt) = p(i, 2);
        pzar(natt) = p(i, 3);
        pear(natt) = p(i, 4);
        xmar(natt) = p(i, 5);
        //Cbz11/11/98end
        //C
        statement_360:;
      }
    }
    //C                ********Fragment the q-qbar jets systems *****
    //C
    jtp(1) = ihnt2(1);
    jtp(2) = ihnt2(3);
    FEM_DO_SAFE(ntp, 1, 2) {
      FEM_DO_SAFE(jjtp, 1, jtp(ntp)) {
        hijfrg(cmn, jjtp, ntp, ierror);
        if (mstu(24) != 0 || ierror > 0) {
          mstu(24) = 0;
          mstu(28) = 0;
          if (ihpr2(10) != 0) {
            //C                  call lulist(2)
            write(6, star), "error occured P&T, repeat the event";
            write(6, star), ntp, jjtp;
            //Clin-6/2009 when this happens, the event will be repeated,
            //C     and another record for the same event number will be written into
            //C     zpc.dat, zpc.res, minijet-initial-beforePropagation.dat,
            //C     parton-initial-afterPropagation.dat, parton-after-coalescence.dat,
            //C     and parton-collisionsHistory.dat.
          }
          goto statement_50;
        }
        //C                        ********check errors
        //C
        nsbst = 1;
        idstr = 92;
        if (ihpr2(21) == 0) {
          luedit(2);
        }
        else {
          statement_381:
          nsbst++;
          if (k(nsbst, 2) < 91 || k(nsbst, 2) > 93) {
            goto statement_381;
          }
          idstr = k(nsbst, 2);
          nsbst++;
        }
        if (frame == "LAB") {
          hboost(cmn);
        }
        //C                ******** boost back to lab frame(if it was in)
        //C
        nftp = nfp(jjtp, 5);
        if (ntp == 2) {
          nftp = 10 + nft(jjtp, 5);
        }
        nsbstr = 0;
        FEM_DO_SAFE(i, nsbst, n) {
          if (k(i, 2) == idstr) {
            nsbstr++;
            goto statement_390;
          }
          k(i, 4) = nsbstr;
          natt++;
          katt(natt, 1) = k(i, 2);
          katt(natt, 2) = nftp;
          katt(natt, 4) = k(i, 1);
          //C                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
          //Clin-4/2008:
          //C                   IF(K(I,3).EQ.0 .OR.
          //C     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
          //C                      KATT(NATT,3)=0
          if (k(i, 3) == 0) {
            katt(natt, 3) = 0;
          }
          else if (k(i, 3) != 0 && k(k(i, 3), 2) == idstr) {
            katt(natt, 3) = 0;
            //Clin-4/2008-end
          }
          else {
            katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i, 3), 4);
          }
          //C       ****** identify the mother particle
          patt(natt, 1) = p(i, 1);
          patt(natt, 2) = p(i, 2);
          patt(natt, 3) = p(i, 3);
          patt(natt, 4) = p(i, 4);
          eatt += p(i, 4);
          //Cbz11/11/98
          //Cbz1/25/99
          //C                   IF (NTP .EQ. 1) THEN
          //C                      GXAR(NATT) = YP(1, jjtp)
          //C                   ELSE
          //C                      GXAR(NATT) = YT(1, jjtp)
          //C                   END IF
          //C                   IF (NTP .EQ. 1) THEN
          //C                      GYAR(NATT) = YP(2, jjtp)
          //C                   ELSE
          //C                      GYAR(NATT) = YT(2, jjtp)
          //C                   END IF
          if (ntp == 1) {
            lsg = jjtp;
          }
          else {
            lsg = jjtp + nsp;
          }
          gxar(natt) = fem::sngl(zt1(lsg));
          gyar(natt) = fem::sngl(zt2(lsg));
          gzar(natt) = fem::sngl(zt3(lsg));
          ftar(natt) = fem::sngl(ataui(lsg));
          //Cbz1/25/99end
          itypar(natt) = k(i, 2);
          pxar(natt) = p(i, 1);
          pyar(natt) = p(i, 2);
          pzar(natt) = p(i, 3);
          pear(natt) = p(i, 4);
          xmar(natt) = p(i, 5);
          //Cbz11/11/98end
          //C
          statement_390:;
        }
      }
    }
    //C     ********Fragment the q-qq related string systems
  }
  //C
  FEM_DO_SAFE(i, 1, ndr) {
    natt++;
    katt(natt, 1) = kfdr(i);
    katt(natt, 2) = 40;
    katt(natt, 3) = 0;
    patt(natt, 1) = pdr(i, 1);
    patt(natt, 2) = pdr(i, 2);
    patt(natt, 3) = pdr(i, 3);
    patt(natt, 4) = pdr(i, 4);
    eatt += pdr(i, 4);
    //Clin-11/11/03     set direct photons positions and time at formation:
    gxar(natt) = rtdr(i, 1);
    gyar(natt) = rtdr(i, 2);
    gzar(natt) = 0.f;
    ftar(natt) = 0.f;
    itypar(natt) = katt(natt, 1);
    pxar(natt) = patt(natt, 1);
    pyar(natt) = patt(natt, 2);
    pzar(natt) = patt(natt, 3);
    pear(natt) = patt(natt, 4);
    xmar(natt) = pdr(i, 5);
  }
  //C
  //C                        ********store the direct-produced particles
  //C
  //Clin-4/19/01 soft3:
  statement_565:
  //C
  dengy = eatt / (ihnt2(1) * hint1(6) + ihnt2(3) * hint1(7)) - 1.0f;
  if (fem::abs(dengy) > hipr1(43) && ihpr2(20) != 0 && ihpr2(21) == 0) {
    if (ihpr2(10) != 0) {
      write(6, star), "Energy not conserved, repeat the event";
    }
    //C                call lulist(1)
    write(6, star), "violated:EATT(GeV),NATT,B(fm)=", eatt, natt, bimp;
    goto statement_50;
  }
  write(6, star), "satisfied:EATT(GeV),NATT,B(fm)=", eatt, natt, bimp;
  write(6, star), " ";
  //C
  //Clin-4/2012 write out initial transverse positions of initial nucleons:
  write(94, star), iaevt, miss, ihnt2(1), ihnt2(3), bimp;
  FEM_DO_SAFE(jp, 1, ihnt2(1)) {
    //Clin-12/2012 write out present and original flavor code of nucleons:
    //C           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP),
    //C     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
    write(94, format_243), yp(1, jp) + 0.5f * bb * fem::cos(phirp),
      yp(2, jp) + 0.5f * bb * fem::sin(phirp), jp, nfp(jp, 5), yp(3,
      jp), nfp(jp, 3), nfp(jp, 4);
  }
  FEM_DO_SAFE(jt, 1, ihnt2(3)) {
    //C target nucleon # has a minus sign for distinction from projectile:
    //Clin-12/2012 write out present and original flavor code of nucleons:
    //C           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP),
    //C     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
    write(94, format_243), yt(1, jt) - 0.5f * bb * fem::cos(phirp),
      yt(2, jt) - 0.5f * bb * fem::sin(phirp), -jt, nft(jt, 5), yt(3,
      jt), nft(jt, 3), nft(jt, 4);
  }
  //Clin-12/2012 write out present and original flavor code of nucleons:
  //C 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
  //Clin-4/2012-end
  //C
}

typedef float (*fnkick_function_pointer)(common&, float const&);

float
fnkick(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  //
  //Cc      SAVE /HPARNT/
  return_value = 1.0f / (x + fem::pow2(hipr1(19))) / (x + fem::pow2(
    hipr1(20))) / (1 + fem::exp((fem::sqrt(x) - hipr1(20)) / 0.4f));
  return return_value;
}

typedef float (*fnkc2_function_pointer)(common&, float const&);

float
fnkc2(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  //
  //Cc      SAVE /HPARNT/
  return_value = x * fem::exp(-2.0f * x / hipr1(42));
  return return_value;
}

typedef float (*fnstru_function_pointer)(common&, float const&);

float
fnstru(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  //Cc      SAVE /HPARNT/
  return_value = fem::pow((1.0f - x), hipr1(44)) / fem::pow((
    fem::pow2(x) + fem::pow2(hipr1(45)) / fem::pow2(hint1(1))), hipr1(
    46));
  return return_value;
}

typedef float (*fnstrm_function_pointer)(common&, float const&);

float
fnstrm(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  //Cc      SAVE /HPARNT/
  return_value = 1.0f / fem::pow((fem::pow2((1.0f - x)) + fem::pow2(
    hipr1(45)) / fem::pow2(hint1(1))), hipr1(46)) / fem::pow((
    fem::pow2(x) + fem::pow2(hipr1(45)) / fem::pow2(hint1(1))), hipr1(
    46));
  return return_value;
}

typedef float (*fnstrs_function_pointer)(common&, float const&);

float
fnstrs(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  //Cc      SAVE /HPARNT/
  return_value = fem::pow((1.0f - x), hipr1(47)) / fem::pow((
    fem::pow2(x) + fem::pow2(hipr1(45)) / fem::pow2(hint1(1))), hipr1(
    48));
  return return_value;
}

float
wdsax(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  // COMMON wood
  float& r = cmn.r;
  float& w = cmn.w;
  //
  //C                             ********THREE PARAMETER WOOD SAXON
  //Cc      SAVE /WOOD/
  return_value = cmn.fnorm * (1.f + w * fem::pow2((x / r))) / (1 +
    fem::exp((x - r) / cmn.d));
  if (w < 0.f) {
    if (x >= r / fem::sqrt(fem::abs(w))) {
      return_value = 0.f;
    }
  }
  return return_value;
}

typedef float (*rwdsax_function_pointer)(common&, float const&);

float
rwdsax(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  return_value = x * x * wdsax(cmn, x);
  return return_value;
}

struct ftot_save
{
  float omg;

  ftot_save() :
    omg(fem::float0)
  {}
};

typedef float (*ftot_function_pointer)(common&, float const&);

float
ftot(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(ftot);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  float& omg = sve.omg;
  //
  //Cc      SAVE /HPARNT/
  omg = omg0(cmn, x) * (hipr1(30) + hint1(11)) / hipr1(31) / 2.0f;
  return_value = 2.0f * (1.0f - fem::exp(-omg));
  return return_value;
}

struct fhin_save
{
  float omg;

  fhin_save() :
    omg(fem::float0)
  {}
};

typedef float (*fhin_function_pointer)(common&, float const&);

float
fhin(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(fhin);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  float& omg = sve.omg;
  //
  //Cc      SAVE /HPARNT/
  omg = omg0(cmn, x) * (hipr1(30) + hint1(11)) / hipr1(31) / 2.0f;
  return_value = 1.0f - fem::exp(-2.0f * omg);
  return return_value;
}

struct ftotjt_save
{
  float omg;

  ftotjt_save() :
    omg(fem::float0)
  {}
};

typedef float (*ftotjt_function_pointer)(common&, float const&);

float
ftotjt(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(ftotjt);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  float& omg = sve.omg;
  //
  //Cc      SAVE /HPARNT/
  omg = omg0(cmn, x) * hint1(11) / hipr1(31) / 2.0f;
  return_value = 1.0f - fem::exp(-2.0f * omg);
  return return_value;
}

struct ftotrg_save
{
  float omg;

  ftotrg_save() :
    omg(fem::float0)
  {}
};

typedef float (*ftotrg_function_pointer)(common&, float const&);

float
ftotrg(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(ftotrg);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  float& omg = sve.omg;
  //
  //Cc      SAVE /HPARNT/
  omg = omg0(cmn, x) * hint1(60) / hipr1(31) / 2.0f;
  return_value = 1.0f - fem::exp(-2.0f * omg);
  return return_value;
}

struct sgmin_save
{
  float ga;
  int i;
  float z;

  sgmin_save() :
    ga(fem::float0),
    i(fem::int0),
    z(fem::float0)
  {}
};

float
sgmin(
  common& cmn,
  int const& n)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(sgmin);
  float& ga = sve.ga;
  int& i = sve.i;
  float& z = sve.z;
  ga = 0.f;
  if (n <= 2) {
    goto statement_20;
  }
  FEM_DO_SAFE(i, 1, n - 1) {
    z = i;
    ga += fem::alog(z);
  }
  statement_20:
  return_value = ga;
  return return_value;
}

struct fnjet_save
{
  float c0;
  float omg1;

  fnjet_save() :
    c0(fem::float0),
    omg1(fem::float0)
  {}
};

typedef float (*fnjet_function_pointer)(common&, float const&);

float
fnjet(
  common& cmn,
  float const& x)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(fnjet);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  // COMMON njet
  int& n = static_cast<common_njet&>(cmn).n;
  //
  // SAVE
  float& c0 = sve.c0;
  float& omg1 = sve.omg1;
  //
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /NJET/
  omg1 = omg0(cmn, x) * hint1(11) / hipr1(31);
  //Clin-8/2015 could cause IEEE_UNDERFLOW, does not seem to affect results:
  c0 = fem::exp(n * fem::alog(omg1) - sgmin(cmn, n + 1));
  if (n == 0) {
    c0 = 1.0f - fem::exp(-2.0f * omg0(cmn, x) * hipr1(30) / hipr1(31) / 2.0f);
  }
  return_value = c0 * fem::exp(-omg1);
  return return_value;
}

struct gauss1_save
{
  float aa;
  float bb;
  float c1;
  float c2;
  float identifier_const;
  float delta;
  int i;
  float s16;
  float s8;
  float u;
  arr<float> w;
  arr<float> x;
  float y;

  gauss1_save() :
    aa(fem::float0),
    bb(fem::float0),
    c1(fem::float0),
    c2(fem::float0),
    identifier_const(fem::float0),
    delta(fem::float0),
    i(fem::int0),
    s16(fem::float0),
    s8(fem::float0),
    u(fem::float0),
    w(dimension(12), fem::fill0),
    x(dimension(12), fem::fill0),
    y(fem::float0)
  {}
};

//C
//C*********GAUSSIAN ONE-DIMENSIONAL INTEGRATION PROGRAM*************
//C
float
gauss1(
  common& cmn,
  fhin_function_pointer f,
  float const& a,
  float const& b,
  float const& eps)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(gauss1);
  common_write write(cmn);
  float& aa = sve.aa;
  float& bb = sve.bb;
  float& c1 = sve.c1;
  float& c2 = sve.c2;
  float& identifier_const = sve.identifier_const;
  float& delta = sve.delta;
  int& i = sve.i;
  float& s16 = sve.s16;
  float& s8 = sve.s8;
  float& u = sve.u;
  arr_ref<float> w(sve.w, dimension(12));
  arr_ref<float> x(sve.x, dimension(12));
  float& y = sve.y;
  if (is_called_first_time) {
    identifier_const = 1.0e-12f;
    {
      static const float values[] = {
        0.1012285f, .2223810f, .3137067f, .3623838f, .0271525f,
          .0622535f, 0.0951585f, .1246290f, .1495960f, .1691565f,
          .1826034f, .1894506f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        w;
    }
    {
      static const float values[] = {
        0.9602899f, .7966665f, .5255324f, .1834346f, .9894009f,
          .9445750f, 0.8656312f, .7554044f, .6178762f, .4580168f,
          .2816036f, .0950125f
      };
      fem::data_of_type<float>(FEM_VALUES_AND_SIZE),
        x;
    }
  }
  //C
  delta = identifier_const * fem::abs(a - b);
  return_value = 0.0f;
  aa = a;
  statement_5:
  y = b - aa;
  if (fem::abs(y) <= delta) {
    return return_value;
  }
  statement_2:
  bb = aa + y;
  c1 = 0.5f * (aa + bb);
  c2 = c1 - aa;
  s8 = 0.0f;
  s16 = 0.0f;
  FEM_DO_SAFE(i, 1, 4) {
    u = x(i) * c2;
    s8 += w(i) * (f(cmn, c1 + u) + f(cmn, c1 - u));
  }
  FEM_DO_SAFE(i, 5, 12) {
    u = x(i) * c2;
    s16 += w(i) * (f(cmn, c1 + u) + f(cmn, c1 - u));
  }
  s8 = s8 * c2;
  s16 = s16 * c2;
  if (fem::abs(s16 - s8) > eps * (1.f + fem::abs(s16))) {
    goto statement_4;
  }
  return_value += s16;
  aa = bb;
  goto statement_5;
  statement_4:
  y = 0.5f * y;
  if (fem::abs(y) > delta) {
    goto statement_2;
  }
  write(6, "(1x,'GAUSS1....TOO HIGH ACURACY REQUIRED')");
  return_value = 0.0f;
  return return_value;
}

struct hifun_save
{
  float fnorm;
  int j;
  float xdd;

  hifun_save() :
    fnorm(fem::float0),
    j(fem::int0),
    xdd(fem::float0)
  {}
};

//C
//C The next three subroutines are for Monte Carlo generation
//C according to a given function FHB. One calls first HIFUN
//C with assigned channel number I, low and up limits. Then to
//C generate the distribution one can call HIRND(I) which gives
//C you a random number generated according to the given function.
//C
void
hifun(
  common& cmn,
  int const& i,
  float const& xmin,
  float const& xmax,
  fnkc2_function_pointer fhb)
{
  FEM_CMN_SVE(hifun);
  // COMMON hijhb
  arr_ref<float, 2> rr(cmn.rr, dimension(10, 201));
  arr_ref<float, 2> xx(cmn.xx, dimension(10, 201));
  //
  // SAVE
  float& fnorm = sve.fnorm;
  int& j = sve.j;
  float& xdd = sve.xdd;
  //
  //Cc      SAVE /HIJHB/
  fnorm = gauss1(cmn, fhb, xmin, xmax, 0.001f);
  FEM_DO_SAFE(j, 1, 201) {
    xx(i, j) = xmin + (xmax - xmin) * (j - 1) / 200.0f;
    xdd = xx(i, j);
    rr(i, j) = gauss1(cmn, fhb, xmin, xdd, 0.001f) / fnorm;
  }
}

struct hijwds_save
{
  float a;
  arr<float> dd;
  float fgaus;
  int i;
  arr<int> iaa;
  arr<float> rr;
  arr<float> ww;
  float xlow;

  hijwds_save() :
    a(fem::float0),
    dd(dimension(20), fem::fill0),
    fgaus(fem::float0),
    i(fem::int0),
    iaa(dimension(20), fem::fill0),
    rr(dimension(20), fem::fill0),
    ww(dimension(20), fem::fill0),
    xlow(fem::float0)
  {}
};

//C
//C ********************************************************
//C ************************              WOOD-SAX
void
hijwds(
  common& cmn,
  int const& ia,
  int const& idh,
  float& xhigh)
{
  FEM_CMN_SVE(hijwds);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  // COMMON wood
  float& r = cmn.r;
  float& d = cmn.d;
  float& fnorm = cmn.fnorm;
  float& w = cmn.w;
  //
  // SAVE
  float& a = sve.a;
  arr_ref<float> dd(sve.dd, dimension(20));
  float& fgaus = sve.fgaus;
  int& i = sve.i;
  arr_ref<int> iaa(sve.iaa, dimension(20));
  arr_ref<float> rr(sve.rr, dimension(20));
  arr_ref<float> ww(sve.ww, dimension(20));
  float& xlow = sve.xlow;
  //
  if (is_called_first_time) {
    {
      fem::data_values data;
      data.values, 2, 4, 12, 16, 27, 32, 40, 56;
      data.values, 63, 93, 184, 197, 208, 7*datum(0.f);
      data, iaa;
    }
    {
      fem::data_values data;
      data.values, 0.01f, .964f, 2.355f, 2.608f, 2.84f, 3.458f, 3.766f, 3.971f;
      data.values, 4.214f, 4.87f, 6.51f, 6.38f, 6.624f, 7*datum(0.f);
      data, rr;
    }
    {
      fem::data_values data;
      data.values, 0.5882f, .322f, .522f, .513f, .569f, .61f, .586f, .5935f;
      data.values, .586f, .573f, .535f, .535f, .549f, 7*datum(0.f);
      data, dd;
    }
    fem::data((values, 0.0f, .517f, -0.149f, -0.051f, 0.f, -0.208f,
      -0.161f, 13*datum(0.f))), ww;
  }
  //C     SETS UP HISTOGRAM IDH WITH RADII FOR
  //C     NUCLEUS IA DISTRIBUTED ACCORDING TO THREE PARAM WOOD SAXON
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /WOOD/
  //C        DIMENSION IAA(20),RR(20),DD(20),WW(20),RMS(20)
  //C
  //C   PARAMETERS OF SPECIAL NUCLEI FROM ATOMIC DATA AND NUC DATA TABLES
  //C     VOL 14, 5-6 1974
  //C        DATA RMS/2.11,1.71,2.46,2.73,3.05,3.247,3.482,3.737,3.925,4.31,
  //C     1        5.42,5.33,5.521,7*0./
  //C
  a = ia;
  //C
  //C                 ********SET WOOD-SAX PARAMS FIRST  AS IN DATE ET AL
  d = 0.54f;
  //C                        ********D IS WOOD SAX DIFFUSE PARAM IN FM
  r = 1.19f * fem::pow(a, (1.f / 3.f)) - 1.61f * fem::pow(a, (-1.f / 3.f));
  //C                         ********R IS RADIUS PARAM
  w = 0.f;
  //C                 ********W IS The third of three WOOD-SAX PARAM
  //C
  //C                      ********CHECK TABLE FOR SPECIAL CASES
  FEM_DO_SAFE(i, 1, 13) {
    if (ia == iaa(i)) {
      r = rr(i);
      d = dd(i);
      w = ww(i);
      //Clin RS not used                              RS=RMS(I)
    }
  }
  //C                             ********FNORM is the normalize factor
  fnorm = 1.0f;
  xlow = 0.f;
  xhigh = r + 12.f * d;
  if (w <  - 0.01f) {
    if (xhigh > r / fem::sqrt(fem::abs(w))) {
      xhigh = r / fem::sqrt(fem::abs(w));
    }
  }
  fgaus = gauss1(cmn, rwdsax, xlow, xhigh, 0.001f);
  fnorm = 1.f / fgaus;
  //C
  if (idh == 1) {
    hint1(72) = r;
    hint1(73) = d;
    hint1(74) = w;
    hint1(75) = fnorm / 4.0f / hipr1(40);
  }
  else if (idh == 2) {
    hint1(76) = r;
    hint1(77) = d;
    hint1(78) = w;
    hint1(79) = fnorm / 4.0f / hipr1(40);
  }
  //C
  //C             NOW SET UP HBOOK FUNCTIONS IDH FOR  R**2*RHO(R)
  //C             THESE HISTOGRAMS ARE USED TO GENERATE RANDOM RADII
  hifun(cmn, idh, xlow, xhigh, rwdsax);
}

double
subcr1(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = 4.e0 / 9.e0 * (1.e0 + fem::pow2(u)) / fem::pow2(t);
  return return_value;
}

double
subcr2(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = 4.e0 / 9.e0 * (fem::pow2(t) + fem::pow2(u));
  return return_value;
}

double
subcr3(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = 4.e0 / 9.e0 * (fem::pow2(t) + fem::pow2(u) + (1.e0 +
    fem::pow2(u)) / fem::pow2(t) - 2.e0 * fem::pow2(u) / 3.e0 / t);
  return return_value;
}

double
subcr4(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = 8.e0 / 3.e0 * (fem::pow2(t) + fem::pow2(u)) * (
    4.e0 / 9.e0 / t / u - 1.e0);
  return return_value;
}

double
subcr5(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = 3.e0 / 8.e0 * (fem::pow2(t) + fem::pow2(u)) * (
    4.e0 / 9.e0 / t / u - 1.e0);
  return return_value;
}

double
subcr6(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = (1.e0 + fem::pow2(u)) * (1.e0 / fem::pow2(t) -
    4.e0 / u / 9.e0);
  return return_value;
}

double
subcr7(
  double const& t,
  double const& u)
{
  double return_value = fem::double0;
  return_value = 9.e0 / 2.e0 * (3.e0 - t * u - u / fem::pow2(t) - t /
    fem::pow2(u));
  return return_value;
}

struct gmre_save
{
  double z;

  gmre_save() :
    z(fem::double0)
  {}
};

double
gmre(
  common& cmn,
  double const& x)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(gmre);
  double& z = sve.z;
  z = x;
  if (x > 3.0e0) {
    goto statement_10;
  }
  z = x + 3.e0;
  statement_10:
  return_value = 0.5e0 * fem::dlog(2.e0 * 3.14159265e0 / z) + z *
    fem::dlog(z) - z + fem::dlog(1.e0 + 1.e0 / 12.e0 / z + 1.e0 /
    288.e0 / fem::pow2(z) - 139.e0 / 51840.e0 / fem::pow3(z) - 571.e0 /
    2488320.e0 / fem::pow4(z));
  if (z == x) {
    goto statement_20;
  }
  return_value = return_value - fem::dlog(z - 1.e0) - fem::dlog(z -
    2.e0) - fem::dlog(z - 3.e0);
  statement_20:
  return return_value;
}

struct parton_save
{
  double aax;
  double ag;
  double aphg;
  double aphs;
  double as;
  double at1;
  double at2;
  double at3;
  double at4;
  double b12;
  double b34;
  double bg;
  double bs;
  double btag;
  double btas;
  double cag;
  double cas;
  double cnd;
  double cnud;
  double dlam;
  double fs1;
  double fs2;
  double fud1;
  double fud2;
  double gmd;
  double gmg;
  double gms;
  double gmud;
  int i;
  double q0;
  double rrx;
  double s;

  parton_save() :
    aax(fem::double0),
    ag(fem::double0),
    aphg(fem::double0),
    aphs(fem::double0),
    as(fem::double0),
    at1(fem::double0),
    at2(fem::double0),
    at3(fem::double0),
    at4(fem::double0),
    b12(fem::double0),
    b34(fem::double0),
    bg(fem::double0),
    bs(fem::double0),
    btag(fem::double0),
    btas(fem::double0),
    cag(fem::double0),
    cas(fem::double0),
    cnd(fem::double0),
    cnud(fem::double0),
    dlam(fem::double0),
    fs1(fem::double0),
    fs2(fem::double0),
    fud1(fem::double0),
    fud2(fem::double0),
    gmd(fem::double0),
    gmg(fem::double0),
    gms(fem::double0),
    gmud(fem::double0),
    i(fem::int0),
    q0(fem::double0),
    rrx(fem::double0),
    s(fem::double0)
  {}
};

void
parton(
  common& cmn,
  arr_ref<double, 2> f,
  double const& x1,
  double const& x2,
  double const& qq)
{
  FEM_CMN_SVE(parton);
  f(dimension(2, 7));
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  int& ipcrs = cmn.ipcrs;
  double& dshadow = cmn.dshadow;
  int& ishadow = cmn.ishadow;
  //
  double& aax = sve.aax;
  double& ag = sve.ag;
  double& aphg = sve.aphg;
  double& aphs = sve.aphs;
  double& as = sve.as;
  double& at1 = sve.at1;
  double& at2 = sve.at2;
  double& at3 = sve.at3;
  double& at4 = sve.at4;
  double& b12 = sve.b12;
  double& b34 = sve.b34;
  double& bg = sve.bg;
  double& bs = sve.bs;
  double& btag = sve.btag;
  double& btas = sve.btas;
  double& cag = sve.cag;
  double& cas = sve.cas;
  double& cnd = sve.cnd;
  double& cnud = sve.cnud;
  double& dlam = sve.dlam;
  double& fs1 = sve.fs1;
  double& fs2 = sve.fs2;
  double& fud1 = sve.fud1;
  double& fud2 = sve.fud2;
  double& gmd = sve.gmd;
  double& gmg = sve.gmg;
  double& gms = sve.gms;
  double& gmud = sve.gmud;
  int& i = sve.i;
  double& q0 = sve.q0;
  double& rrx = sve.rrx;
  double& s = sve.s;
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /NJET/
  //Clin-7/2009:
  dlam = fem::dble(hipr1(15));
  q0 = fem::dble(hipr1(16));
  s = fem::dlog(fem::dlog(qq / fem::pow2(dlam)) / fem::dlog(fem::pow2(
    q0) / fem::pow2(dlam)));
  if (ihpr2(7) == 2) {
    goto statement_200;
  }
  //C*******************************************************
  at1 = 0.419e0 + 0.004e0 * s - 0.007e0 * fem::pow2(s);
  at2 = 3.460e0 + 0.724e0 * s - 0.066e0 * fem::pow2(s);
  gmud = 4.40e0 - 4.86e0 * s + 1.33e0 * fem::pow2(s);
  at3 = 0.763e0 - 0.237e0 * s + 0.026e0 * fem::pow2(s);
  at4 = 4.00e0 + 0.627e0 * s - 0.019e0 * fem::pow2(s);
  gmd = -0.421e0 * s + 0.033e0 * fem::pow2(s);
  //C*******************************************************
  cas = 1.265e0 - 1.132e0 * s + 0.293e0 * fem::pow2(s);
  as = -0.372e0 * s - 0.029e0 * fem::pow2(s);
  bs = 8.05e0 + 1.59e0 * s - 0.153e0 * fem::pow2(s);
  aphs = 6.31e0 * s - 0.273e0 * fem::pow2(s);
  btas = -10.5e0 * s - 3.17e0 * fem::pow2(s);
  gms = 14.7e0 * s + 9.80e0 * fem::pow2(s);
  //C********************************************************
  //C        CAC=0.135*S-0.075*S**2
  //C        AC=-0.036-0.222*S-0.058*S**2
  //C        BC=6.35+3.26*S-0.909*S**2
  //C        APHC=-3.03*S+1.50*S**2
  //C        BTAC=17.4*S-11.3*S**2
  //C        GMC=-17.9*S+15.6*S**2
  //C***********************************************************
  cag = 1.56e0 - 1.71e0 * s + 0.638e0 * fem::pow2(s);
  ag = -0.949e0 * s + 0.325e0 * fem::pow2(s);
  bg = 6.0e0 + 1.44e0 * s - 1.05e0 * fem::pow2(s);
  aphg = 9.0e0 - 7.19e0 * s + 0.255e0 * fem::pow2(s);
  btag = -16.5e0 * s + 10.9e0 * fem::pow2(s);
  gmg = 15.3e0 * s - 10.1e0 * fem::pow2(s);
  goto statement_300;
  //C********************************************************
  statement_200:
  at1 = 0.374e0 + 0.014e0 * s;
  at2 = 3.33e0 + 0.753e0 * s - 0.076e0 * fem::pow2(s);
  gmud = 6.03e0 - 6.22e0 * s + 1.56e0 * fem::pow2(s);
  at3 = 0.761e0 - 0.232e0 * s + 0.023e0 * fem::pow2(s);
  at4 = 3.83e0 + 0.627e0 * s - 0.019e0 * fem::pow2(s);
  gmd = -0.418e0 * s + 0.036e0 * fem::pow2(s);
  //C************************************
  cas = 1.67e0 - 1.92e0 * s + 0.582e0 * fem::pow2(s);
  as = -0.273e0 * s - 0.164e0 * fem::pow2(s);
  bs = 9.15e0 + 0.530e0 * s - 0.763e0 * fem::pow2(s);
  aphs = 15.7e0 * s - 2.83e0 * fem::pow2(s);
  btas = -101.0e0 * s + 44.7e0 * fem::pow2(s);
  gms = 223.0e0 * s - 117.0e0 * fem::pow2(s);
  //C*********************************
  //C        CAC=0.067*S-0.031*S**2
  //C        AC=-0.120-0.233*S-0.023*S**2
  //C        BC=3.51+3.66*S-0.453*S**2
  //C        APHC=-0.474*S+0.358*S**2
  //C        BTAC=9.50*S-5.43*S**2
  //C        GMC=-16.6*S+15.5*S**2
  //C**********************************
  cag = 0.879e0 - 0.971e0 * s + 0.434e0 * fem::pow2(s);
  ag = -1.16e0 * s + 0.476e0 * fem::pow2(s);
  bg = 4.0e0 + 1.23e0 * s - 0.254e0 * fem::pow2(s);
  aphg = 9.0e0 - 5.64e0 * s - 0.817e0 * fem::pow2(s);
  btag = -7.54e0 * s + 5.50e0 * fem::pow2(s);
  gmg = -0.596e0 * s + 1.26e0 * fem::pow2(s);
  //C*********************************
  statement_300:
  b12 = fem::dexp(gmre(cmn, at1) + gmre(cmn, at2 + 1.e0) - gmre(cmn,
    at1 + at2 + 1.e0));
  b34 = fem::dexp(gmre(cmn, at3) + gmre(cmn, at4 + 1.e0) - gmre(cmn,
    at3 + at4 + 1.e0));
  cnud = 3.e0 / b12 / (1.e0 + gmud * at1 / (at1 + at2 + 1.e0));
  cnd = 1.e0 / b34 / (1.e0 + gmd * at3 / (at3 + at4 + 1.e0));
  //C********************************************************
  //C        FUD=X*(U+D)
  //C        FS=X*2(UBAR+DBAR+SBAR)  AND UBAR=DBAR=SBAR
  //C*******************************************************
  fud1 = cnud * fem::pow(x1, at1) * fem::pow((1.e0 - x1), at2) * (
    1.e0 + gmud * x1);
  fs1 = cas * fem::pow(x1, as) * fem::pow((1.e0 - x1), bs) * (1.e0 +
    aphs * x1 + btas * fem::pow2(x1) + gms * fem::pow3(x1));
  f(1, 3) = cnd * fem::pow(x1, at3) * fem::pow((1.e0 - x1), at4) * (
    1.e0 + gmd * x1) + fs1 / 6.e0;
  f(1, 1) = fud1 - f(1, 3) + fs1 / 3.e0;
  f(1, 2) = fs1 / 6.e0;
  f(1, 4) = fs1 / 6.e0;
  f(1, 5) = fs1 / 6.e0;
  f(1, 6) = fs1 / 6.e0;
  f(1, 7) = cag * fem::pow(x1, ag) * fem::pow((1.e0 - x1), bg) * (
    1.e0 + aphg * x1 + btag * fem::pow2(x1) + gmg * fem::pow3(x1));
  //C
  fud2 = cnud * fem::pow(x2, at1) * fem::pow((1.e0 - x2), at2) * (
    1.e0 + gmud * x2);
  fs2 = cas * fem::pow(x2, as) * fem::pow((1.e0 - x2), bs) * (1.e0 +
    aphs * x2 + btas * fem::pow2(x2) + gms * fem::pow3(x2));
  f(2, 3) = cnd * fem::pow(x2, at3) * fem::pow((1.e0 - x2), at4) * (
    1.e0 + gmd * x2) + fs2 / 6.e0;
  f(2, 1) = fud2 - f(2, 3) + fs2 / 3.e0;
  f(2, 2) = fs2 / 6.e0;
  f(2, 4) = fs2 / 6.e0;
  f(2, 5) = fs2 / 6.e0;
  f(2, 6) = fs2 / 6.e0;
  f(2, 7) = cag * fem::pow(x2, ag) * fem::pow((1.e0 - x2), bg) * (
    1.e0 + aphg * x2 + btag * fem::pow2(x2) + gmg * fem::pow3(x2));
  //C***********Nuclear effect on the structure function****************
  //C
  if (ihpr2(6) == 1 && ihnt2(1) > 1) {
    aax = 1.193e0 * fem::dble(fem::pow(fem::alog(fem::ffloat(ihnt2(1))),
      0.16666666f));
    rrx = aax * (fem::pow3(x1) - 1.2e0 * fem::pow2(x1) + 0.21e0 *
      x1) + 1.e0 + fem::dble(1.079f * (fem::pow(fem::ffloat(ihnt2(1)),
      0.33333333f) - 1.0f)) / fem::dble(fem::alog(fem::ffloat(ihnt2(
      1)) + 1.0f)) * fem::dsqrt(x1) * fem::dexp(-fem::pow2(x1) /
      0.01e0);
    //Clin-8/2015 DEXP() above may cause IEEE_UNDERFLOW,
    //C     does not seem to affect results.
    //C     &          /DLOG(IHNT2(1)+1.0D0)*(DSQRT(X1)*DEXP(-X1**2/0.01)
    //Clin-7/2009 enable users to modify nuclear shadowing:
    if (ishadow == 1) {
      rrx = 1.e0 + dshadow * (rrx - 1.e0);
    }
    if (ipcrs == 1 || ipcrs == 3) {
      rrx = fem::dexp(-fem::pow2(x1) / 0.01e0);
    }
    //Clin-7/2009:
    if ((ipcrs == 1 || ipcrs == 3) && ishadow == 1) {
      rrx = fem::dexp(-fem::pow2(x1) / 0.01e0) * dshadow;
    }
    FEM_DO_SAFE(i, 1, 7) {
      f(1, i) = rrx * f(1, i);
    }
  }
  if (ihpr2(6) == 1 && ihnt2(3) > 1) {
    aax = 1.193e0 * fem::dble(fem::pow(fem::alog(fem::ffloat(ihnt2(3))),
      0.16666666f));
    rrx = aax * (fem::pow3(x2) - 1.2e0 * fem::pow2(x2) + 0.21e0 *
      x2) + 1.e0 + fem::dble(1.079f * (fem::pow(fem::ffloat(ihnt2(3)),
      0.33333f) - 1.0f)) / fem::dble(fem::alog(fem::ffloat(ihnt2(
      3)) + 1.0f)) * fem::dsqrt(x2) * fem::dexp(-fem::pow2(x2) /
      0.01e0);
    //C     &         /DLOG(IHNT2(3)+1.0D0)*DSQRT(X2)*DEXP(-X2**2/0.01)
    //Clin-7/2009:
    if (ishadow == 1) {
      rrx = 1.e0 + dshadow * (rrx - 1.e0);
    }
    if (ipcrs == 2 || ipcrs == 3) {
      rrx = fem::dexp(-fem::pow2(x2) / 0.01e0);
    }
    //Clin-7/2009:
    if ((ipcrs == 2 || ipcrs == 3) && ishadow == 1) {
      rrx = fem::dexp(-fem::pow2(x2) / 0.01e0) * dshadow;
    }
    FEM_DO_SAFE(i, 1, 7) {
      f(2, i) = rrx * f(2, i);
    }
  }
  //C
}

struct g_save
{
  double af;
  double aph;
  double dlam;
  arr<double, 2> f;
  double g11;
  double g12;
  double g13;
  double g2;
  double g31;
  double g32;
  double g4;
  double g5;
  double g61;
  double g62;
  double g7;
  double ss;
  double t;
  double u;
  double x1;
  double x2;
  double xt;
  double z;

  g_save() :
    af(fem::double0),
    aph(fem::double0),
    dlam(fem::double0),
    f(dimension(2, 7), fem::fill0),
    g11(fem::double0),
    g12(fem::double0),
    g13(fem::double0),
    g2(fem::double0),
    g31(fem::double0),
    g32(fem::double0),
    g4(fem::double0),
    g5(fem::double0),
    g61(fem::double0),
    g62(fem::double0),
    g7(fem::double0),
    ss(fem::double0),
    t(fem::double0),
    u(fem::double0),
    x1(fem::double0),
    x2(fem::double0),
    xt(fem::double0),
    z(fem::double0)
  {}
};

double
g(
  common& cmn,
  double const& y1,
  double const& y2,
  double const& pt2)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(g);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  double& af = sve.af;
  double& aph = sve.aph;
  double& dlam = sve.dlam;
  arr_ref<double, 2> f(sve.f, dimension(2, 7));
  double& g11 = sve.g11;
  double& g12 = sve.g12;
  double& g13 = sve.g13;
  double& g2 = sve.g2;
  double& g31 = sve.g31;
  double& g32 = sve.g32;
  double& g4 = sve.g4;
  double& g5 = sve.g5;
  double& g61 = sve.g61;
  double& g62 = sve.g62;
  double& g7 = sve.g7;
  double& ss = sve.ss;
  double& t = sve.t;
  double& u = sve.u;
  double& x1 = sve.x1;
  double& x2 = sve.x2;
  double& xt = sve.xt;
  double& z = sve.z;
  //
  //Cc      SAVE /HPARNT/
  xt = 2.e0 * fem::dsqrt(pt2) / fem::dble(hint1(1));
  x1 = 0.5e0 * xt * (fem::dexp(y1) + fem::dexp(y2));
  x2 = 0.5e0 * xt * (fem::dexp(-y1) + fem::dexp(-y2));
  z = fem::dsqrt(1.e0 - fem::pow2(xt) / x1 / x2);
  ss = x1 * x2 * fem::pow2(fem::dble(hint1(1)));
  t = -(1.e0 - z) / 2.e0;
  u = -(1.e0 + z) / 2.e0;
  af = 3.e0;
  dlam = fem::dble(hipr1(15));
  aph = 12.e0 * 3.1415926e0 / (33.e0 - 2.e0 * af) / fem::dlog(pt2 /
    fem::pow2(dlam));
  //C
  parton(cmn, f, x1, x2, pt2);
  //C
  g11 = ((f(1, 1) + f(1, 2)) * (f(2, 3) + f(2, 4) + f(2, 5) + f(2,
    6)) + (f(1, 3) + f(1, 4)) * (f(2, 5) + f(2, 6))) * subcr1(t, u);
  //C
  g12 = ((f(2, 1) + f(2, 2)) * (f(1, 3) + f(1, 4) + f(1, 5) + f(1,
    6)) + (f(2, 3) + f(2, 4)) * (f(1, 5) + f(1, 6))) * subcr1(u, t);
  //C
  g13 = (f(1, 1) * f(2, 1) + f(1, 2) * f(2, 2) + f(1, 3) * f(2, 3) + f(1,
    4) * f(2, 4) + f(1, 5) * f(2, 5) + f(1, 6) * f(2, 6)) * (subcr1(u,
    t) + subcr1(t, u) - 8.e0 / t / u / 27.e0);
  //C
  g2 = (af - 1) * (f(1, 1) * f(2, 2) + f(2, 1) * f(1, 2) + f(1, 3) * f(2,
    4) + f(2, 3) * f(1, 4) + f(1, 5) * f(2, 6) + f(2, 5) * f(1, 6)) * subcr2(t,
    u);
  //C
  g31 = (f(1, 1) * f(2, 2) + f(1, 3) * f(2, 4) + f(1, 5) * f(2, 6)) * subcr3(t,
    u);
  g32 = (f(2, 1) * f(1, 2) + f(2, 3) * f(1, 4) + f(2, 5) * f(1, 6)) * subcr3(u,
    t);
  //C
  g4 = (f(1, 1) * f(2, 2) + f(2, 1) * f(1, 2) + f(1, 3) * f(2, 4) + f(2,
    3) * f(1, 4) + f(1, 5) * f(2, 6) + f(2, 5) * f(1, 6)) * subcr4(t,
    u);
  //C
  g5 = af * f(1, 7) * f(2, 7) * subcr5(t, u);
  //C
  g61 = f(1, 7) * (f(2, 1) + f(2, 2) + f(2, 3) + f(2, 4) + f(2, 5) + f(2,
    6)) * subcr6(t, u);
  g62 = f(2, 7) * (f(1, 1) + f(1, 2) + f(1, 3) + f(1, 4) + f(1, 5) + f(1,
    6)) * subcr6(u, t);
  //C
  g7 = f(1, 7) * f(2, 7) * subcr7(t, u);
  //C
  return_value = (g11 + g12 + g13 + g2 + g31 + g32 + g4 + g5 + g61 +
    g62 + g7) * fem::dble(hipr1(17)) * 3.14159e0 * fem::pow2(aph) /
    fem::pow2(ss);
  return return_value;
}

struct fjet_save
{
  double pt2;
  double xt;
  double y1;
  double y2;
  double ymn2;
  double ymx1;
  double ymx2;

  fjet_save() :
    pt2(fem::double0),
    xt(fem::double0),
    y1(fem::double0),
    y2(fem::double0),
    ymn2(fem::double0),
    ymx1(fem::double0),
    ymx2(fem::double0)
  {}
};

typedef double (*fjet_function_pointer)(common&, arr_cref<double>,
  double const&);

double
fjet(
  common& cmn,
  arr_cref<double> x,
  double const& /* wgt */)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(fjet);
  x(dimension(10));
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  double& pt2 = sve.pt2;
  double& xt = sve.xt;
  double& y1 = sve.y1;
  double& y2 = sve.y2;
  double& ymn2 = sve.ymn2;
  double& ymx1 = sve.ymx1;
  double& ymx2 = sve.ymx2;
  //
  //Cc      SAVE /HPARNT/
  pt2 = fem::dble(fem::pow2(hint1(1)) / 4.0f - fem::pow2(hipr1(8))) *
    x(1) + fem::pow2(fem::dble(hipr1(8)));
  xt = 2.0e0 * fem::dsqrt(pt2) / fem::dble(hint1(1));
  ymx1 = fem::dlog(1.0e0 / xt + fem::dsqrt(1.0e0 / fem::pow2(xt) - 1.0e0));
  y1 = 2.0e0 * ymx1 * x(2) - ymx1;
  ymx2 = fem::dlog(2.0e0 / xt - fem::dexp(y1));
  ymn2 = fem::dlog(2.0e0 / xt - fem::dexp(-y1));
  y2 = (ymx2 + ymn2) * x(3) - ymn2;
  return_value = 2.0e0 * ymx1 * (ymx2 + ymn2) * fem::dble(fem::pow2(
    hint1(1)) / 4.0f - fem::pow2(hipr1(8))) * g(cmn, y1, y2, pt2) /
    2.0e0;
  return return_value;
}

struct ghvq_save
{
  double af;
  double aph;
  double dlam;
  arr<double, 2> f;
  double ggg;
  double gqq;
  double ss;
  double x1;
  double x2;
  double xt;

  ghvq_save() :
    af(fem::double0),
    aph(fem::double0),
    dlam(fem::double0),
    f(dimension(2, 7), fem::fill0),
    ggg(fem::double0),
    gqq(fem::double0),
    ss(fem::double0),
    x1(fem::double0),
    x2(fem::double0),
    xt(fem::double0)
  {}
};

double
ghvq(
  common& cmn,
  double const& y1,
  double const& y2,
  double const& amt2)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(ghvq);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  double& af = sve.af;
  double& aph = sve.aph;
  double& dlam = sve.dlam;
  arr_ref<double, 2> f(sve.f, dimension(2, 7));
  double& ggg = sve.ggg;
  double& gqq = sve.gqq;
  double& ss = sve.ss;
  double& x1 = sve.x1;
  double& x2 = sve.x2;
  double& xt = sve.xt;
  //
  //Cc      SAVE /HPARNT/
  xt = 2.0e0 * fem::dsqrt(amt2) / fem::dble(hint1(1));
  x1 = 0.5e0 * xt * (fem::dexp(y1) + fem::dexp(y2));
  x2 = 0.5e0 * xt * (fem::dexp(-y1) + fem::dexp(-y2));
  ss = x1 * x2 * fem::pow2(fem::dble(hint1(1)));
  af = 4.0e0;
  if (ihpr2(18) != 0) {
    af = 5.0e0;
  }
  dlam = fem::dble(hipr1(15));
  aph = 12.0e0 * 3.1415926e0 / (33.e0 - 2.e0 * af) / fem::dlog(amt2 /
    fem::pow2(dlam));
  //C
  parton(cmn, f, x1, x2, amt2);
  //C
  gqq = 4.e0 * (dcosh(y1 - y2) + fem::pow2(fem::dble(hipr1(7))) / amt2) / (
    1.e0 + dcosh(y1 - y2)) / 9.e0 * (f(1, 1) * f(2, 2) + f(1, 2) * f(2,
    1) + f(1, 3) * f(2, 4) + f(1, 4) * f(2, 3) + f(1, 5) * f(2, 6) + f(1,
    6) * f(2, 5));
  ggg = (8.e0 * dcosh(y1 - y2) - 1.e0) * (dcosh(y1 - y2) + 2.e0 *
    fem::pow2(fem::dble(hipr1(7))) / amt2 - 2.e0 * fem::pow4(fem::dble(
    hipr1(7))) / fem::pow2(amt2)) / (1.e0 + dcosh(y1 - y2)) / 24.e0 * f(1,
    7) * f(2, 7);
  //C
  return_value = (gqq + ggg) * fem::dble(hipr1(23)) * 3.14159e0 *
    fem::pow2(aph) / fem::pow2(ss);
  return return_value;
}

struct gphotn_save
{
  double af;
  double aph;
  double aphem;
  double dlam;
  arr<double, 2> f;
  double g11;
  double g12;
  double g2;
  double ss;
  double t;
  double u;
  double x1;
  double x2;
  double xt;
  double z;

  gphotn_save() :
    af(fem::double0),
    aph(fem::double0),
    aphem(fem::double0),
    dlam(fem::double0),
    f(dimension(2, 7), fem::fill0),
    g11(fem::double0),
    g12(fem::double0),
    g2(fem::double0),
    ss(fem::double0),
    t(fem::double0),
    u(fem::double0),
    x1(fem::double0),
    x2(fem::double0),
    xt(fem::double0),
    z(fem::double0)
  {}
};

double
gphotn(
  common& cmn,
  double const& y1,
  double const& y2,
  double const& pt2)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(gphotn);
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  double& af = sve.af;
  double& aph = sve.aph;
  double& aphem = sve.aphem;
  double& dlam = sve.dlam;
  arr_ref<double, 2> f(sve.f, dimension(2, 7));
  double& g11 = sve.g11;
  double& g12 = sve.g12;
  double& g2 = sve.g2;
  double& ss = sve.ss;
  double& t = sve.t;
  double& u = sve.u;
  double& x1 = sve.x1;
  double& x2 = sve.x2;
  double& xt = sve.xt;
  double& z = sve.z;
  //
  //Cc      SAVE /HPARNT/
  xt = 2.e0 * fem::dsqrt(pt2) / fem::dble(hint1(1));
  x1 = 0.5e0 * xt * (fem::dexp(y1) + fem::dexp(y2));
  x2 = 0.5e0 * xt * (fem::dexp(-y1) + fem::dexp(-y2));
  z = fem::dsqrt(1.e0 - fem::pow2(xt) / x1 / x2);
  ss = x1 * x2 * fem::pow2(fem::dble(hint1(1)));
  t = -(1.e0 - z) / 2.e0;
  u = -(1.e0 + z) / 2.e0;
  af = 3.e0;
  dlam = fem::dble(hipr1(15));
  aph = 12.e0 * 3.1415926e0 / (33.e0 - 2.e0 * af) / fem::dlog(pt2 /
    fem::pow2(dlam));
  aphem = 1.e0 / 137.e0;
  //C
  parton(cmn, f, x1, x2, pt2);
  //C
  g11 = -(fem::pow2(u) + 1.e0) / u / 3.e0 * f(1, 7) * (4.e0 * f(2,
    1) + 4.e0 * f(2, 2) + f(2, 3) + f(2, 4) + f(2, 5) + f(2, 6)) /
    9.e0;
  g12 = -(fem::pow2(t) + 1.e0) / t / 3.e0 * f(2, 7) * (4.e0 * f(1,
    1) + 4.e0 * f(1, 2) + f(1, 3) + f(1, 4) + f(1, 5) + f(1, 6)) /
    9.e0;
  g2 = 8.e0 * (fem::pow2(u) + fem::pow2(t)) / u / t / 9.e0 * (4.e0 * f(1,
    1) * f(2, 2) + 4.e0 * f(1, 2) * f(2, 1) + f(1, 3) * f(2, 4) + f(1,
    4) * f(2, 3) + f(1, 5) * f(2, 6) + f(1, 6) * f(2, 5)) / 9.e0;
  //C
  return_value = (g11 + g12 + g2) * fem::dble(hipr1(23)) *
    3.14159e0 * aph * aphem / fem::pow2(ss);
  return return_value;
}

struct fjetrg_save
{
  double am2;
  double amt2;
  double gtrig;
  double pt2;
  float ptmax;
  float ptmin;
  double xt;
  double y1;
  double y2;
  double ymn2;
  double ymx1;
  double ymx2;

  fjetrg_save() :
    am2(fem::double0),
    amt2(fem::double0),
    gtrig(fem::double0),
    pt2(fem::double0),
    ptmax(fem::float0),
    ptmin(fem::float0),
    xt(fem::double0),
    y1(fem::double0),
    y2(fem::double0),
    ymn2(fem::double0),
    ymx1(fem::double0),
    ymx2(fem::double0)
  {}
};

typedef double (*fjetrg_function_pointer)(common&, arr_cref<double>,
  double const&);

double
fjetrg(
  common& cmn,
  arr_cref<double> x,
  double const& /* wgt */)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(fjetrg);
  x(dimension(10));
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  //
  // SAVE
  double& am2 = sve.am2;
  double& amt2 = sve.amt2;
  double& gtrig = sve.gtrig;
  double& pt2 = sve.pt2;
  float& ptmax = sve.ptmax;
  float& ptmin = sve.ptmin;
  double& xt = sve.xt;
  double& y1 = sve.y1;
  double& y2 = sve.y2;
  double& ymn2 = sve.ymn2;
  double& ymx1 = sve.ymx1;
  double& ymx2 = sve.ymx2;
  //
  //Cc      SAVE /HPARNT/
  ptmin = fem::abs(hipr1(10)) - 0.25f;
  ptmin = fem::max(ptmin, hipr1(8));
  am2 = 0.e0;
  if (ihpr2(3) == 3) {
    am2 = fem::dble(fem::pow2(hipr1(7)));
    ptmin = fem::max(0.0f, hipr1(10));
  }
  ptmax = fem::abs(hipr1(10)) + 0.25f;
  if (hipr1(10) <= 0.0f) {
    ptmax = hint1(1) / 2.0f - fem::sngl(am2);
  }
  if (ptmax <= ptmin) {
    ptmax = ptmin + 0.25f;
  }
  pt2 = fem::dble(fem::pow2(ptmax) - fem::pow2(ptmin)) * x(1) +
    fem::pow2(fem::dble(ptmin));
  amt2 = pt2 + am2;
  xt = 2.0e0 * fem::dsqrt(amt2) / fem::dble(hint1(1));
  ymx1 = fem::dlog(1.0e0 / xt + fem::dsqrt(1.0e0 / fem::pow2(xt) - 1.0e0));
  y1 = 2.0e0 * ymx1 * x(2) - ymx1;
  ymx2 = fem::dlog(2.0e0 / xt - fem::dexp(y1));
  ymn2 = fem::dlog(2.0e0 / xt - fem::dexp(-y1));
  y2 = (ymx2 + ymn2) * x(3) - ymn2;
  if (ihpr2(3) == 3) {
    gtrig = 2.0e0 * ghvq(cmn, y1, y2, amt2);
  }
  else if (ihpr2(3) == 2) {
    gtrig = 2.0e0 * gphotn(cmn, y1, y2, pt2);
  }
  else {
    gtrig = g(cmn, y1, y2, pt2);
  }
  return_value = 2.0e0 * ymx1 * (ymx2 + ymn2) * fem::dble(fem::pow2(
    ptmax) - fem::pow2(ptmin)) * gtrig / 2.0e0;
  return return_value;
}

struct aran9_save
{
  int i;

  aran9_save() :
    i(fem::int0)
  {}
};

void
aran9(
  common& cmn,
  arr_ref<float> qran,
  int const& ndim)
{
  FEM_CMN_SVE(aran9);
  qran(dimension(10));
  // SAVE
  int& i = sve.i;
  //
  FEM_DO_SAFE(i, 1, ndim) {
    qran(i) = ranart(cmn.num1);
  }
}

struct vegas_save
{
  double alph;
  double calls;
  arr<double, 2> d;
  arr<double, 2> di;
  double dr;
  arr<double> dt;
  double dv2g;
  arr<double> dx;
  double dxg;
  double f2;
  double f2b;
  double fb;
  int i;
  arr<int> ia;
  int j;
  int k;
  arr<int> kg;
  int mds;
  int nd;
  int ndm;
  int ndmx;
  int ng;
  int npg;
  double one;
  arr<float> qran;
  arr<double> r;
  double rc;
  double ti2;
  double wgt;
  arr<double> x;
  arr<double> xin;
  double xjac;
  double xn;
  double xnd;
  double xo;

  vegas_save() :
    alph(fem::double0),
    calls(fem::double0),
    d(dimension(50, 10), fem::fill0),
    di(dimension(50, 10), fem::fill0),
    dr(fem::double0),
    dt(dimension(10), fem::fill0),
    dv2g(fem::double0),
    dx(dimension(10), fem::fill0),
    dxg(fem::double0),
    f2(fem::double0),
    f2b(fem::double0),
    fb(fem::double0),
    i(fem::int0),
    ia(dimension(10), fem::fill0),
    j(fem::int0),
    k(fem::int0),
    kg(dimension(10), fem::fill0),
    mds(fem::int0),
    nd(fem::int0),
    ndm(fem::int0),
    ndmx(fem::int0),
    ng(fem::int0),
    npg(fem::int0),
    one(fem::double0),
    qran(dimension(10), fem::fill0),
    r(dimension(50), fem::fill0),
    rc(fem::double0),
    ti2(fem::double0),
    wgt(fem::double0),
    x(dimension(10), fem::fill0),
    xin(dimension(50), fem::fill0),
    xjac(fem::double0),
    xn(fem::double0),
    xnd(fem::double0),
    xo(fem::double0)
  {}
};

//C*******************************************************************
//C
//C*******************************************************************
//C   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
//C      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
//C*******************************************************************
//C
void
vegas(
  common& cmn,
  fjet_function_pointer fxn,
  double& avgi,
  double& sd,
  double& chi2a)
{
  FEM_CMN_SVE(vegas);
  common_write write(cmn);
  arr_cref<double> xl(cmn.xl, dimension(10));
  arr_cref<double> xu(cmn.xu, dimension(10));
  int& ndim = cmn.ndim;
  int& ncall = cmn.ncall;
  arr_ref<double, 2> xi(cmn.xi, dimension(50, 10));
  double& si = cmn.si;
  double& si2 = cmn.si2;
  double& swgt = cmn.swgt;
  double& schi = cmn.schi;
  int& ndo = cmn.ndo;
  int& it = cmn.it;
  double& f = cmn.f;
  double& ti = cmn.ti;
  double& tsi = cmn.tsi;
  //
  double& alph = sve.alph;
  double& calls = sve.calls;
  arr_ref<double, 2> d(sve.d, dimension(50, 10));
  arr_ref<double, 2> di(sve.di, dimension(50, 10));
  double& dr = sve.dr;
  arr_ref<double> dt(sve.dt, dimension(10));
  double& dv2g = sve.dv2g;
  arr_ref<double> dx(sve.dx, dimension(10));
  double& dxg = sve.dxg;
  double& f2 = sve.f2;
  double& f2b = sve.f2b;
  double& fb = sve.fb;
  int& i = sve.i;
  arr_ref<int> ia(sve.ia, dimension(10));
  int& j = sve.j;
  int& k = sve.k;
  arr_ref<int> kg(sve.kg, dimension(10));
  int& mds = sve.mds;
  int& nd = sve.nd;
  int& ndm = sve.ndm;
  int& ndmx = sve.ndmx;
  int& ng = sve.ng;
  int& npg = sve.npg;
  double& one = sve.one;
  arr_ref<float> qran(sve.qran, dimension(10));
  arr_ref<double> r(sve.r, dimension(50));
  double& rc = sve.rc;
  double& ti2 = sve.ti2;
  double& wgt = sve.wgt;
  arr_ref<double> x(sve.x, dimension(10));
  arr_ref<double> xin(sve.xin, dimension(50));
  double& xjac = sve.xjac;
  double& xn = sve.xn;
  double& xnd = sve.xnd;
  double& xo = sve.xo;
  if (is_called_first_time) {
    ndmx = 50;
    alph = 1.5e0;
    one = 1.e0;
    mds = -1;
  }
  //Cc      SAVE /BVEG1/
  //Cc      SAVE /BVEG2/
  //Cc      SAVE /BVEG3/
  //C      REAL*4 QRAN(10)
  //C
  ndo = 1;
  FEM_DO_SAFE(j, 1, ndim) {
    xi(1, j) = one;
  }
  //C
  // UNHANDLED: ENTRY vegas1(fxn,avgi,sd,chi2a)
  //C         - INITIALIZES CUMMULATIVE VARIABLES, BUT NOT GRID
  it = 0;
  si = 0.e0;
  si2 = si;
  swgt = si;
  schi = si;
  //C
  // UNHANDLED: ENTRY vegas2(fxn,avgi,sd,chi2a)
  //C         - NO INITIALIZATION
  nd = ndmx;
  ng = 1;
  if (mds == 0) {
    goto statement_2;
  }
  ng = fem::fint(fem::pow((fem::real(ncall) / 2.f), (1.f / fem::real(ndim))));
  mds = 1;
  if ((2 * ng - ndmx) < 0) {
    goto statement_2;
  }
  mds = -1;
  npg = ng / ndmx + 1;
  nd = ng / npg;
  ng = npg * nd;
  statement_2:
  k = fem::pow(ng, ndim);
  npg = ncall / k;
  if (npg < 2) {
    npg = 2;
  }
  calls = npg * k;
  dxg = one / ng;
  dv2g = fem::pow2((calls * fem::pow(dxg, ndim))) / npg / npg / (npg - one);
  xnd = nd;
  ndm = nd - 1;
  dxg = dxg * xnd;
  xjac = one / calls;
  FEM_DO_SAFE(j, 1, ndim) {
    //C***this is the line 50
    dx(j) = xu(j) - xl(j);
    xjac = xjac * dx(j);
  }
  //C
  //C   REBIN PRESERVING BIN DENSITY
  //C
  if (nd == ndo) {
    goto statement_8;
  }
  rc = ndo / xnd;
  FEM_DO_SAFE(j, 1, ndim) {
    k = 0;
    xn = 0.e0;
    dr = xn;
    i = k;
    statement_4:
    k++;
    dr += one;
    xo = xn;
    xn = xi(k, j);
    statement_5:
    if (rc > dr) {
      goto statement_4;
    }
    i++;
    dr = dr - rc;
    xin(i) = xn - (xn - xo) * dr;
    if (i < ndm) {
      goto statement_5;
    }
    FEM_DO_SAFE(i, 1, ndm) {
      xi(i, j) = xin(i);
    }
    xi(nd, j) = one;
  }
  ndo = nd;
  //C
  statement_8:
  //C      IF(NPRN.NE.0) WRITE(16,200) NDIM,CALLS,IT,ITMX,ACC,MDS,ND
  //C     1                           ,(XL(J),XU(J),J=1,NDIM)
  //C
  // UNHANDLED: ENTRY vegas3(fxn,avgi,sd,chi2a)
  //C         - MAIN INTEGRATION LOOP
  statement_9:
  it++;
  ti = 0.e0;
  tsi = ti;
  FEM_DO_SAFE(j, 1, ndim) {
    kg(j) = 1;
    FEM_DO_SAFE(i, 1, nd) {
      d(i, j) = ti;
      di(i, j) = ti;
    }
  }
  //C
  statement_11:
  fb = 0.e0;
  f2b = fb;
  k = 0;
  statement_12:
  k++;
  aran9(cmn, qran, ndim);
  wgt = xjac;
  FEM_DO_SAFE(j, 1, ndim) {
    xn = fem::dble(fem::ffloat(kg(j)) - qran(j)) * dxg + one;
    //C*****this is the line 100
    ia(j) = fem::fint(xn);
    if (ia(j) > 1) {
      goto statement_13;
    }
    xo = xi(ia(j), j);
    rc = (xn - ia(j)) * xo;
    goto statement_14;
    statement_13:
    xo = xi(ia(j), j) - xi(ia(j) - 1, j);
    rc = xi(ia(j) - 1, j) + (xn - ia(j)) * xo;
    statement_14:
    x(j) = xl(j) + rc * dx(j);
    wgt = wgt * xo * xnd;
  }
  //C
  f = wgt;
  f = f * fxn(cmn, x, wgt);
  f2 = f * f;
  fb += f;
  f2b += f2;
  FEM_DO_SAFE(j, 1, ndim) {
    di(ia(j), j) += f;
    if (mds >= 0) {
      d(ia(j), j) += f2;
    }
  }
  if (k < npg) {
    goto statement_12;
  }
  //C
  f2b = fem::dsqrt(f2b * npg);
  f2b = (f2b - fb) * (f2b + fb);
  ti += fb;
  tsi += f2b;
  if (mds >= 0) {
    goto statement_18;
  }
  FEM_DO_SAFE(j, 1, ndim) {
    d(ia(j), j) += f2b;
  }
  statement_18:
  k = ndim;
  statement_19:
  kg(k) = fem::mod(kg(k), ng) + 1;
  if (kg(k) != 1) {
    goto statement_11;
  }
  k = k - 1;
  if (k > 0) {
    goto statement_19;
  }
  //C
  //C   FINAL RESULTS FOR THIS ITERATION
  //C
  tsi = tsi * dv2g;
  ti2 = ti * ti;
  wgt = ti2 / (tsi + 1.0e-37);
  si += ti * wgt;
  si2 += ti2;
  swgt += wgt;
  swgt += 1.0e-37;
  si2 += 1.0e-37;
  schi += ti2 * wgt;
  avgi = si / swgt;
  sd = swgt * it / si2;
  chi2a = sd * (schi / swgt - avgi * avgi) / fem::dble(fem::ffloat(it) - .999f);
  sd = fem::dsqrt(one / sd);
  //C****this is the line 150
  if (cmn.nprn == 0) {
    goto statement_21;
  }
  tsi = fem::dsqrt(tsi);
  //C      WRITE(16,201) IT,TI,TSI,AVGI,SD,CHI2A
  //C      IF(NPRN.GE.0) GO TO 21
  //C      DO 20 J=1,NDIM
  //C20    WRITE(16,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
  //C
  //C   REFINE GRID
  //C
  statement_21:
  FEM_DO_SAFE(j, 1, ndim) {
    xo = d(1, j);
    xn = d(2, j);
    d(1, j) = (xo + xn) / 2.e0;
    dt(j) = d(1, j);
    FEM_DO_SAFE(i, 2, ndm) {
      d(i, j) = xo + xn;
      xo = xn;
      xn = d(i + 1, j);
      d(i, j) = (d(i, j) + xn) / 3.e0;
      dt(j) += d(i, j);
    }
    d(nd, j) = (xn + xo) / 2.e0;
    dt(j) += d(nd, j);
  }
  //C
  FEM_DO_SAFE(j, 1, ndim) {
    rc = 0.e0;
    FEM_DO_SAFE(i, 1, nd) {
      r(i) = 0.e0;
      if (dt(j) >= 1.0e18) {
        write(6, star), "************** A SINGULARITY >1.0D18";
        //C      WRITE(5,1111)
        //C1111  FORMAT(1X,'**************IMPORTANT NOTICE***************')
        //C      WRITE(5,1112)
        //C1112  FORMAT(1X,'THE INTEGRAND GIVES RISE A SINGULARITY >1.0D18')
        //C      WRITE(5,1113)
        //C1113  FORMAT(1X,'PLEASE CHECK THE INTEGRAND AND THE LIMITS')
        //C      WRITE(5,1114)
        //C1114  FORMAT(1X,'**************END NOTICE*************')
      }
      if (d(i, j) <= 1.0e-18) {
        goto statement_24;
      }
      xo = dt(j) / d(i, j);
      r(i) = fem::pow(((xo - one) / xo / fem::dlog(xo)), alph);
      statement_24:
      rc += r(i);
    }
    rc = rc / xnd;
    k = 0;
    xn = 0.e0;
    dr = xn;
    i = k;
    statement_25:
    k++;
    dr += r(k);
    xo = xn;
    //C****this is the line 200
    xn = xi(k, j);
    statement_26:
    if (rc > dr) {
      goto statement_25;
    }
    i++;
    dr = dr - rc;
    xin(i) = xn - (xn - xo) * dr / (r(k) + 1.0e-30);
    if (i < ndm) {
      goto statement_26;
    }
    FEM_DO_SAFE(i, 1, ndm) {
      xi(i, j) = xin(i);
    }
    xi(nd, j) = one;
  }
  //C
  if (it < cmn.itmx && cmn.acc * fem::dabs(avgi) < sd) {
    goto statement_9;
  }
  //C200   FORMAT('0INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,'  NCALL=',F8.0
  //C     1    /28X,'  IT=',I5,'  ITMX=',I5/28X,'  ACC=',G9.3
  //C     2    /28X,'  MDS=',I3,'   ND=',I4/28X,'  (XL,XU)=',
  //C     3    (T40,'( ',G12.6,' , ',G12.6,' )'))
  //C201   FORMAT(///' INTEGRATION BY VEGAS' / '0ITERATION NO.',I3,
  //C     1    ':   INTEGRAL =',G14.8/21X,'STD DEV  =',G10.4 /
  //C     2    ' ACCUMULATED RESULTS:   INTEGRAL =',G14.8 /
  //C     3    24X,'STD DEV  =',G10.4 / 24X,'CHI**2 PER IT''N =',G10.4)
  //C202   FORMAT('0DATA FOR AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
  //C     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
  //C     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
  //C     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
}

struct crsjet_save
{
  double avgi;
  double chi2a;
  double sd;

  crsjet_save() :
    avgi(fem::double0),
    chi2a(fem::double0),
    sd(fem::double0)
  {}
};

//C
//C        THIS PROGRAM IS TO CALCULATE THE JET CROSS SECTION
//C        THE INTEGRATION IS DONE BY USING VEGAS
//C
void
crsjet(
  common& cmn)
{
  FEM_CMN_SVE(crsjet);
  // COMMON hparnt
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  // COMMON njet
  int& ipcrs = cmn.ipcrs;
  //
  // SAVE
  double& avgi = sve.avgi;
  double& chi2a = sve.chi2a;
  double& sd = sve.sd;
  //
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /NJET/
  //Cc      SAVE /BVEG1/
  //Cc      SAVE /BVEG2/
  //Cc      SAVE /BVEG3/
  //Cc      SAVE /SEDVAX/
  //C
  //C************************
  //C        NCALL give the number of inner-iteration, ITMX
  //C       gives the limit of out-iteration. Nprn is an option
  //C       ( 1: print the integration process. 0: do not print)
  //C
  cmn.ndim = 3;
  ipcrs = 0;
  vegas(cmn, fjet, avgi, sd, chi2a);
  hint1(14) = fem::sngl(avgi) / 2.5682f;
  if (ihpr2(6) == 1 && ihnt2(1) > 1) {
    ipcrs = 1;
    vegas(cmn, fjet, avgi, sd, chi2a);
    hint1(15) = fem::sngl(avgi) / 2.5682f;
  }
  if (ihpr2(6) == 1 && ihnt2(3) > 1) {
    ipcrs = 2;
    vegas(cmn, fjet, avgi, sd, chi2a);
    hint1(16) = fem::sngl(avgi) / 2.5682f;
  }
  if (ihpr2(6) == 1 && ihnt2(1) > 1 && ihnt2(3) > 1) {
    ipcrs = 3;
    vegas(cmn, fjet, avgi, sd, chi2a);
    hint1(17) = fem::sngl(avgi) / 2.5682f;
  }
  //C                ********Total inclusive jet cross section(Pt>P0)
  //C
  if (ihpr2(3) != 0) {
    ipcrs = 0;
    vegas(cmn, fjetrg, avgi, sd, chi2a);
    hint1(61) = fem::sngl(avgi) / 2.5682f;
    if (ihpr2(6) == 1 && ihnt2(1) > 1) {
      ipcrs = 1;
      vegas(cmn, fjetrg, avgi, sd, chi2a);
      hint1(62) = fem::sngl(avgi) / 2.5682f;
    }
    if (ihpr2(6) == 1 && ihnt2(3) > 1) {
      ipcrs = 2;
      vegas(cmn, fjetrg, avgi, sd, chi2a);
      hint1(63) = fem::sngl(avgi) / 2.5682f;
    }
    if (ihpr2(6) == 1 && ihnt2(1) > 1 && ihnt2(3) > 1) {
      ipcrs = 3;
      vegas(cmn, fjetrg, avgi, sd, chi2a);
      hint1(64) = fem::sngl(avgi) / 2.5682f;
    }
  }
  //C                        ********cross section of trigger jet
  //C
}

struct hijcrs_save
{
  float aphx1;
  float aphx2;
  int i;

  hijcrs_save() :
    aphx1(fem::float0),
    aphx2(fem::float0),
    i(fem::int0)
  {}
};

void
hijcrs(
  common& cmn)
{
  FEM_CMN_SVE(hijcrs);
  // COMMON hparnt
  arr_ref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  // COMMON njet
  int& n = static_cast<common_njet&>(cmn).n;
  //
  // SAVE
  float& aphx1 = sve.aphx1;
  float& aphx2 = sve.aphx2;
  int& i = sve.i;
  //
  //C        THIS IS TO CALCULATE THE CROSS SECTIONS OF JET PRODUCTION AND
  //C        THE TOTAL INELASTIC CROSS SECTIONS.
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /NJET/
  if (hint1(1) >= 10.0f) {
    crsjet(cmn);
  }
  //C                        ********calculate jet cross section(in mb)
  //C
  //Clin-7/2009 these are related to nuclear shadowing:
  aphx1 = hipr1(6) * (fem::pow(ihnt2(1), 0.3333333f) - 1.0f);
  aphx2 = hipr1(6) * (fem::pow(ihnt2(3), 0.3333333f) - 1.0f);
  hint1(11) = hint1(14) - aphx1 * hint1(15) - aphx2 * hint1(16) +
    aphx1 * aphx2 * hint1(17);
  hint1(10) = gauss1(cmn, ftotjt, 0.0f, 20.0f, 0.01f);
  hint1(12) = gauss1(cmn, fhin, 0.0f, 20.0f, 0.01f);
  hint1(13) = gauss1(cmn, ftot, 0.0f, 20.0f, 0.01f);
  hint1(60) = hint1(61) - aphx1 * hint1(62) - aphx2 * hint1(63) +
    aphx1 * aphx2 * hint1(64);
  hint1(59) = gauss1(cmn, ftotrg, 0.0f, 20.0f, 0.01f);
  if (hint1(59) == 0.0f) {
    hint1(59) = hint1(60);
  }
  if (hint1(1) >= 10.0f) {
    FEM_DO_SAFE(i, 0, 20) {
      n = i;
      hint1(80 + i) = gauss1(cmn, fnjet, 0.0f, 20.0f, 0.01f) / hint1(12);
    }
  }
  hint1(10) = hint1(10) * hipr1(31);
  hint1(12) = hint1(12) * hipr1(31);
  hint1(13) = hint1(13) * hipr1(31);
  hint1(59) = hint1(59) * hipr1(31);
  //C                ********Total and Inel cross section are calculated
  //C                        by Gaussian integration.
  if (ihpr2(13) != 0) {
    hipr1(33) = 1.36f * (1.0f + 36.0f / fem::pow2(hint1(1))) *
      fem::alog(0.6f + 0.1f * fem::pow2(hint1(1)));
    hipr1(33) = hipr1(33) / hint1(12);
  }
  //C                ********Parametrized cross section for single
  //C                        diffractive reaction(Goulianos)
}

void
title(
  common& cmn)
{
  common_write write(cmn);
  //C
  //Cc      SAVE /RNDF77/
  //C
  write(6,
    "(/,/,10x,'**************************************************',/,10x,"
    "'*     |      |       _______      /  ------/     *',/,10x,"
    "'*   ----- ------     |_____|     /_/     /       *',/,10x,"
    "'*    |||    /        |_____|      /    / |       *',/,10x,"
    "'*    /| |  /_/       /_______    /_  /    |      *',/,10x,"
    "'*   / |     / /     /  /  / |        -------     *',/,10x,"
    "'*     |    / /|       /  /  |     /     |        *',/,10x,"
    "'*     |   / /  |     /  /  _|    /   -------     *',/,10x,"
    "'*                                                *',/,10x,"
    "'**************************************************',/,10x,"
    "'                      HIJING                      ',/,10x,"
    "'       Heavy Ion Jet INteraction Generator        ',/,10x,"
    "'                        by                        ',/,10x,"
    "'            X. N. Wang  and  M. Gyulassy           ',/,10x,"
    "'             Lawrence Berkeley Laboratory           ',/,/)");
  //Clin-8/15/02 f77:
  //C200        FORMAT(//10X,
  //C     &        '**************************************************'/10X,
  //C     &  '*     |      \       _______      /  ------/     *'/10X,
  //C     &        '*   ----- ------     |_____|     /_/     /       *'/10X,
  //C     &        '*    ||\    /        |_____|      /    / \       *'/10X,
  //C     &        '*    /| \  /_/       /_______    /_  /    \_     *'/10X,
  //C     &        '*   / |     / /     /  /  / |        -------     *'/10X,
  //C     &        '*     |    / /\       /  /  |     /     |        *'/10X,
  //C     &        '*     |   / /  \     /  / \_|    /   -------     *'/10X,
}

struct hijset_save
{
  double dd1;
  double dd2;
  double dd3;
  double dd4;
  fem::str<4> eframe;
  int i;
  int j;
  float rkp;
  float rmax;

  hijset_save() :
    dd1(fem::double0),
    dd2(fem::double0),
    dd3(fem::double0),
    dd4(fem::double0),
    eframe(fem::char0),
    i(fem::int0),
    j(fem::int0),
    rkp(fem::float0),
    rmax(fem::float0)
  {}
};

void
hijset(
  common& cmn,
  float const& efrm,
  str_cref frame,
  str_cref proj,
  str_cref targ,
  int const& iap,
  int const& izp,
  int const& iat,
  int const& izt)
{
  FEM_CMN_SVE(hijset);
  common_write write(cmn);
  arr_ref<float> hipr1(cmn.hipr1, dimension(100));
  arr_cref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  arr_ref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<float, 2> hidat0(cmn.hidat0, dimension(10, 10));
  arr_ref<float> hidat(cmn.hidat, dimension(10));
  arr_ref<int> mstu(cmn.mstu, dimension(200));
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  arr_ref<float> parj(cmn.parj, dimension(200));
  //
  double& dd1 = sve.dd1;
  double& dd2 = sve.dd2;
  double& dd3 = sve.dd3;
  double& dd4 = sve.dd4;
  fem::str<4>& eframe = sve.eframe;
  int& i = sve.i;
  int& j = sve.j;
  float& rkp = sve.rkp;
  float& rmax = sve.rmax;
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HIJDAT/
  //Cc      SAVE /LUDAT1/
  //C
  title(cmn);
  ihnt2(1) = iap;
  ihnt2(2) = izp;
  ihnt2(3) = iat;
  ihnt2(4) = izt;
  ihnt2(5) = 0;
  ihnt2(6) = 0;
  //C
  hint1(8) = fem::max(ulmass(2112), ulmass(2212));
  hint1(9) = hint1(8);
  //C
  if (proj != "A") {
    if (proj == "P") {
      ihnt2(5) = 2212;
    }
    else if (proj == "PBAR") {
      ihnt2(5) = -2212;
    }
    else if (proj == "PI+") {
      ihnt2(5) = 211;
    }
    else if (proj == "PI-") {
      ihnt2(5) = -211;
    }
    else if (proj == "K+") {
      ihnt2(5) = 321;
    }
    else if (proj == "K-") {
      ihnt2(5) = -321;
    }
    else if (proj == "N") {
      ihnt2(5) = 2112;
    }
    else if (proj == "NBAR") {
      ihnt2(5) = -2112;
    }
    else {
      write(6, star), proj, "wrong or unavailable proj name";
      FEM_STOP(0);
    }
    hint1(8) = ulmass(ihnt2(5));
  }
  if (targ != "A") {
    if (targ == "P") {
      ihnt2(6) = 2212;
    }
    else if (targ == "PBAR") {
      ihnt2(6) = -2212;
    }
    else if (targ == "PI+") {
      ihnt2(6) = 211;
    }
    else if (targ == "PI-") {
      ihnt2(6) = -211;
    }
    else if (targ == "K+") {
      ihnt2(6) = 321;
    }
    else if (targ == "K-") {
      ihnt2(6) = -321;
    }
    else if (targ == "N") {
      ihnt2(6) = 2112;
    }
    else if (targ == "NBAR") {
      ihnt2(6) = -2112;
    }
    else {
      write(6, star), targ, "wrong or unavailable targ name";
      FEM_STOP(0);
    }
    hint1(9) = ulmass(ihnt2(6));
  }
  //C
  //C...Switch off decay of pi0, K0S, Lambda, Sigma+-, Xi0-, Omega-.
  if (ihpr2(12) > 0) {
    lugive("MDCY(C221,1)=0");
    //Clin-11/07/00 no K* decays:
    lugive("MDCY(C313,1)=0");
    lugive("MDCY(C-313,1)=0");
    lugive("MDCY(C323,1)=0");
    lugive("MDCY(C-323,1)=0");
    //Clin-1/04/01 no K0 and K0bar decays so K0L and K0S do not appear,
    //C     this way the K/Kbar difference is accounted for exactly:
    lugive("MDCY(C311,1)=0");
    lugive("MDCY(C-311,1)=0");
    //Clin-11/08/00 no Delta decays:
    lugive("MDCY(C1114,1)=0");
    lugive("MDCY(C2114,1)=0");
    lugive("MDCY(C2214,1)=0");
    lugive("MDCY(C2224,1)=0");
    lugive("MDCY(C-1114,1)=0");
    lugive("MDCY(C-2114,1)=0");
    lugive("MDCY(C-2214,1)=0");
    lugive("MDCY(C-2224,1)=0");
    //Clin-11/07/00-end
    //Cbz12/4/98
    lugive("MDCY(C213,1)=0");
    lugive("MDCY(C-213,1)=0");
    lugive("MDCY(C113,1)=0");
    lugive("MDCY(C223,1)=0");
    lugive("MDCY(C333,1)=0");
    //Cbz12/4/98end
    lugive("MDCY(C111,1)=0");
    lugive("MDCY(C310,1)=0");
    lugive("MDCY(C411,1)=0;MDCY(C-411,1)=0");
    lugive("MDCY(C421,1)=0;MDCY(C-421,1)=0");
    lugive("MDCY(C431,1)=0;MDCY(C-431,1)=0");
    lugive("MDCY(C511,1)=0;MDCY(C-511,1)=0");
    lugive("MDCY(C521,1)=0;MDCY(C-521,1)=0");
    lugive("MDCY(C531,1)=0;MDCY(C-531,1)=0");
    lugive("MDCY(C3122,1)=0;MDCY(C-3122,1)=0");
    lugive("MDCY(C3112,1)=0;MDCY(C-3112,1)=0");
    lugive("MDCY(C3212,1)=0;MDCY(C-3212,1)=0");
    lugive("MDCY(C3222,1)=0;MDCY(C-3222,1)=0");
    lugive("MDCY(C3312,1)=0;MDCY(C-3312,1)=0");
    lugive("MDCY(C3322,1)=0;MDCY(C-3322,1)=0");
    lugive("MDCY(C3334,1)=0;MDCY(C-3334,1)=0");
    //Clin-7/2011-no HQ(charm or bottom) decays in order to get net-HQ conservation:
    lugive("MDCY(C441,1)=0");
    lugive("MDCY(C443,1)=0");
    lugive("MDCY(C413,1)=0;MDCY(C-413,1)=0");
    lugive("MDCY(C423,1)=0;MDCY(C-423,1)=0");
    lugive("MDCY(C433,1)=0;MDCY(C-433,1)=0");
    lugive("MDCY(C4112,1)=0;MDCY(C-4112,1)=0");
    lugive("MDCY(C4114,1)=0;MDCY(C-4114,1)=0");
    lugive("MDCY(C4122,1)=0;MDCY(C-4122,1)=0");
    lugive("MDCY(C4212,1)=0;MDCY(C-4212,1)=0");
    lugive("MDCY(C4214,1)=0;MDCY(C-4214,1)=0");
    lugive("MDCY(C4222,1)=0;MDCY(C-4222,1)=0");
    lugive("MDCY(C4224,1)=0;MDCY(C-4224,1)=0");
    lugive("MDCY(C4132,1)=0;MDCY(C-4132,1)=0");
    lugive("MDCY(C4312,1)=0;MDCY(C-4312,1)=0");
    lugive("MDCY(C4314,1)=0;MDCY(C-4314,1)=0");
    lugive("MDCY(C4232,1)=0;MDCY(C-4232,1)=0");
    lugive("MDCY(C4322,1)=0;MDCY(C-4322,1)=0");
    lugive("MDCY(C4324,1)=0;MDCY(C-4324,1)=0");
    lugive("MDCY(C4332,1)=0;MDCY(C-4332,1)=0");
    lugive("MDCY(C4334,1)=0;MDCY(C-4334,1)=0");
    lugive("MDCY(C551,1)=0");
    lugive("MDCY(C553,1)=0");
    lugive("MDCY(C513,1)=0;MDCY(C-513,1)=0");
    lugive("MDCY(C523,1)=0;MDCY(C-523,1)=0");
    lugive("MDCY(C533,1)=0;MDCY(C-533,1)=0");
    lugive("MDCY(C5112,1)=0;MDCY(C-5112,1)=0");
    lugive("MDCY(C5114,1)=0;MDCY(C-5114,1)=0");
    lugive("MDCY(C5122,1)=0;MDCY(C-5122,1)=0");
    lugive("MDCY(C5212,1)=0;MDCY(C-5212,1)=0");
    lugive("MDCY(C5214,1)=0;MDCY(C-5214,1)=0");
    lugive("MDCY(C5222,1)=0;MDCY(C-5222,1)=0");
    lugive("MDCY(C5224,1)=0;MDCY(C-5224,1)=0");
    //Clin-7/2011-end
  }
  mstu(12) = 0;
  mstu(21) = 1;
  if (ihpr2(10) == 0) {
    mstu(22) = 0;
    mstu(25) = 0;
    mstu(26) = 0;
  }
  //C
  //Clin    parj(41) and (42) are a, b parameters in Lund, read from input.ampt:
  //C        PARJ(41)=HIPR1(3)
  //C        PARJ(42)=HIPR1(4)
  //C        PARJ(41)=2.2
  //C        PARJ(42)=0.5
  //C
  //Clin  2 popcorn parameters read from input.ampt:
  //C        IHPR2(11) = 3
  //C        PARJ(5) = 0.5
  mstj(12) = ihpr2(11);
  //C
  //Clin  parj(21) gives the mean gaussian width for hadron Pt:
  parj(21) = hipr1(2);
  //Clin  parj(2) is gamma_s=P(s)/P(u), kappa propto 1/b/(2+a) assumed.
  rkp = hipr1(4) * (2 + hipr1(3)) / parj(42) / (2 + parj(41));
  parj(2) = fem::pow(parj(2), (1.f / rkp));
  parj(21) = parj(21) * fem::sqrt(rkp);
  //Clin-10/31/00 update when string tension is changed:
  hipr1(2) = parj(21);
  //C
  //Clin-4/2015: set upper limit for gamma_s=P(s)/P(u) to 0.4
  //C     (to limit strangeness enhancement when string tension is strongly
  //C     increased due to using a very low value of parameter b in Lund
  //C     symmetric splitting function as done in arXiv:1403.6321):
  parj(2) = fem::min(parj(2), 0.4f);
  //C
  //C                        ******** set up for jetset
  if (frame == "LAB") {
    dd1 = fem::dble(efrm);
    dd2 = fem::dble(hint1(8));
    dd3 = fem::dble(hint1(9));
    hint1(1) = fem::sqrt(fem::pow2(hint1(8)) + 2.0f * hint1(9) *
      efrm + fem::pow2(hint1(9)));
    dd4 = fem::dsqrt(fem::pow2(dd1) - fem::pow2(dd2)) / (dd1 + dd3);
    hint1(2) = fem::sngl(dd4);
    hint1(3) = 0.5f * fem::sngl(fem::dlog((1.e0 + dd4) / (1.e0 - dd4)));
    dd4 = fem::dsqrt(fem::pow2(dd1) - fem::pow2(dd2)) / dd1;
    hint1(4) = 0.5f * fem::sngl(fem::dlog((1.e0 + dd4) / (1.e0 - dd4)));
    hint1(5) = 0.0f;
    hint1(6) = efrm;
    hint1(7) = hint1(9);
  }
  else if (frame == "CMS") {
    hint1(1) = efrm;
    hint1(2) = 0.0f;
    hint1(3) = 0.0f;
    dd1 = fem::dble(hint1(1));
    dd2 = fem::dble(hint1(8));
    dd3 = fem::dble(hint1(9));
    dd4 = fem::dsqrt(1.e0 - 4.e0 * fem::pow2(dd2) / fem::pow2(dd1));
    hint1(4) = 0.5f * fem::sngl(fem::dlog((1.e0 + dd4) / (1.e0 - dd4)));
    dd4 = fem::dsqrt(1.e0 - 4.e0 * fem::pow2(dd3) / fem::pow2(dd1));
    hint1(5) = -0.5f * fem::sngl(fem::dlog((1.e0 + dd4) / (1.e0 - dd4)));
    hint1(6) = hint1(1) / 2.0f;
    hint1(7) = hint1(1) / 2.0f;
  }
  //C                ********define Lorentz transform to lab frame
  //C
  //C                ********calculate the cross sections involved with
  //C                        nucleon collisions.
  if (ihnt2(1) > 1) {
    hijwds(cmn, ihnt2(1), 1, rmax);
    hipr1(34) = rmax;
    //C                        ********set up Wood-Sax distr for proj.
  }
  if (ihnt2(3) > 1) {
    hijwds(cmn, ihnt2(3), 2, rmax);
    hipr1(35) = rmax;
    //C                        ********set up Wood-Sax distr for  targ.
  }
  //C
  i = 0;
  statement_20:
  i++;
  if (i == 10) {
    goto statement_30;
  }
  if (hidat0(10, i) <= hint1(1)) {
    goto statement_20;
  }
  statement_30:
  if (i == 1) {
    i = 2;
  }
  FEM_DO_SAFE(j, 1, 9) {
    hidat(j) = hidat0(j, i - 1) + (hidat0(j, i) - hidat0(j, i - 1)) *
      (hint1(1) - hidat0(10, i - 1)) / (hidat0(10, i) - hidat0(10,
      i - 1));
  }
  hipr1(31) = hidat(5);
  hipr1(30) = 2.0f * hidat(5);
  //C
  hijcrs(cmn);
  //C
  if (ihpr2(5) != 0) {
    hifun(cmn, 3, 0.0f, 36.0f, fnkick);
    //C                ********booking for generating pt**2 for pt kick
  }
  hifun(cmn, 7, 0.0f, 6.0f, fnkc2);
  hifun(cmn, 4, 0.0f, 1.0f, fnstru);
  hifun(cmn, 5, 0.0f, 1.0f, fnstrm);
  hifun(cmn, 6, 0.0f, 1.0f, fnstrs);
  //C                ********booking for x distribution of valence quarks
  eframe = "Ecm";
  if (frame == "LAB") {
    eframe = "Elab";
  }
  write(6,
    "(10x,'**************************************************',/,10x,'*',48x,"
    "'*',/,10x,'*         HIJING has been initialized at         *',/,10x,'*',"
    "13x,a4,'= ',f10.2,' GeV/n',13x,'*',/,10x,'*',48x,'*',/,10x,'*',8x,'for ',"
    "a4,'(',i3,',',i3,')',' + ',a4,'(',i3,',',i3,')',7x,'*',/,10x,"
    "'**************************************************')"),
    eframe, efrm, proj, ihnt2(1), ihnt2(2), targ, ihnt2(3), ihnt2(4);
}

struct blockdata_hidata_save
{
  int i;
  int j;

  blockdata_hidata_save() :
    i(fem::int0),
    j(fem::int0)
  {}
};

//C
//C***************************************************************
//C
void
blockdata_hidata(
  common& cmn)
{
  FEM_CMN_SVE(blockdata_hidata);
  // COMMON bveg1
  arr_ref<double> xl(cmn.xl, dimension(10));
  arr_ref<double> xu(cmn.xu, dimension(10));
  // COMMON hparnt
  arr_ref<float> hipr1(cmn.hipr1, dimension(100));
  arr_ref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_ref<float> hint1(cmn.hint1, dimension(100));
  arr_ref<int> ihnt2(cmn.ihnt2, dimension(50));
  // COMMON hmain2
  const int maxstr = 150001;
  arr_ref<int, 2> katt(cmn.katt, dimension(maxstr, 4));
  arr_ref<float, 2> patt(cmn.patt, dimension(maxstr, 4));
  // COMMON hstrng
  arr_ref<int, 2> nfp(cmn.nfp, dimension(300, 15));
  arr_ref<float, 2> pp(cmn.pp, dimension(300, 15));
  arr_ref<int, 2> nft(cmn.nft, dimension(300, 15));
  arr_ref<float, 2> pt(cmn.pt, dimension(300, 15));
  // COMMON hjcrdn
  arr_ref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_ref<float, 2> yt(cmn.yt, dimension(3, 300));
  // COMMON hjjet1
  arr_ref<int> npj(cmn.npj, dimension(300));
  arr_ref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_ref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_ref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_ref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_ref<float, 2> pjpe(cmn.pjpe, dimension(300, 500));
  arr_ref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_ref<int> ntj(cmn.ntj, dimension(300));
  arr_ref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_ref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_ref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_ref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_ref<float, 2> pjte(cmn.pjte, dimension(300, 500));
  arr_ref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  // COMMON hjjet2
  arr_ref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_ref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_ref<int, 2> k1sg(cmn.k1sg, dimension(maxstr, 100));
  arr_ref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_ref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_ref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_ref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_ref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_ref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  // COMMON hijdat
  arr_ref<float, 2> hidat0(cmn.hidat0, dimension(10, 10));
  arr_ref<float> hidat(cmn.hidat, dimension(10));
  // COMMON hpint
  arr_ref<float, 2> atco(cmn.atco, dimension(200, 20));
  arr_ref<float> atxs(cmn.atxs, dim1(0, 200));
  //
  if (is_called_first_time) {
    cmn.num1 = 30123984;
    fem::data((values, 10*datum(0.e0))), xl;
    fem::data((values, 10*datum(1.e0))), xu;
    cmn.ncall = 1000;
    cmn.itmx = 100;
    cmn.acc = 0.01f;
    cmn.nprn = 0;
    {
      fem::data_values data;
      data.values, 1.5f, 0.35f, 0.5f, 0.9f, 2.0f, 0.1f, 1.5f, 2.0f;
      data.values, -1.0f, -2.25f, 2.0f, 0.5f, 1.0f, 2.0f, 0.2f, 2.0f;
      data.values, 2.5f, 0.3f, 0.1f, 1.4f, 1.6f, 1.0f, 2.0f, 0.0f;
      data.values, 0.0f, 0.0f, 0.0f, 0.0f, 0.4f, 57.0f, 28.5f, 3.9f;
      data.values, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.14159f;
      data.values, 0.0f, 0.4f, 0.1f, 1.5f, 0.1f, 0.25f, 0.0f, 0.5f;
      data.values, 0.0f, 0.0f, 50*datum(0.0f);
      data, hipr1;
    }
    {
      fem::data_values data;
      data.values, 1, 3, 0, 1, 1, 1, 1, 10;
      data.values, 0, 0, 1, 1, 1, 1, 0, 0;
      data.values, 1, 0, 0, 1, 30*datum(0);
      data, ihpr2;
    }
    fem::data((values, 100*datum(0))), hint1;
    fem::data((values, 50*datum(0))), ihnt2;
    cmn.natt = 0;
    cmn.eatt = 0.0f;
    cmn.jatt = 0;
    cmn.nt = 0;
    cmn.np = 0;
    cmn.n0 = 0;
    cmn.n01 = 0;
    cmn.n10 = 0;
    cmn.n11 = 0;
    fem::data((values, 600004*datum(0))), katt;
    fem::data((values, 600004*datum(0.0f))), patt;
    fem::data((values, 4500*datum(0))), nfp;
    fem::data((values, 4500*datum(0.0f))), pp;
    fem::data((values, 4500*datum(0))), nft;
    fem::data((values, 4500*datum(0.0f))), pt;
    fem::data((values, 900*datum(0.0f))), yp;
    fem::data((values, 900*datum(0.0f))), yt;
    fem::data((values, 300*datum(0))), npj;
    fem::data((values, 150000*datum(0))), kfpj;
    fem::data((values, 150000*datum(0.0f))), pjpx;
    fem::data((values, 150000*datum(0.0f))), pjpy;
    fem::data((values, 150000*datum(0.0f))), pjpz;
    fem::data((values, 150000*datum(0.0f))), pjpe;
    fem::data((values, 150000*datum(0.0f))), pjpm;
    fem::data((values, 300*datum(0))), ntj;
    fem::data((values, 150000*datum(0))), kftj;
    fem::data((values, 150000*datum(0.0f))), pjtx;
    fem::data((values, 150000*datum(0.0f))), pjty;
    fem::data((values, 150000*datum(0.0f))), pjtz;
    fem::data((values, 150000*datum(0.0f))), pjte;
    fem::data((values, 150000*datum(0.0f))), pjtm;
    cmn.nsg = 0;
    fem::data((values, 150001*datum(0))), njsg;
    fem::data((values, 450003*datum(0))), iasg;
    fem::data((values, 15000100*datum(0))), k1sg;
    fem::data((values, 15000100*datum(0))), k2sg;
    fem::data((values, 15000100*datum(0.0f))), pxsg;
    fem::data((values, 15000100*datum(0.0f))), pysg;
    fem::data((values, 15000100*datum(0.0f))), pzsg;
    fem::data((values, 15000100*datum(0.0f))), pesg;
    fem::data((values, 15000100*datum(0.0f))), pmsg;
    cmn.mint4 = 0;
    cmn.mint5 = 0;
    fem::data((values, 4000*datum(0.0f))), atco;
    fem::data((values, 201*datum(0.0f))), atxs;
    {
      static const float values[] = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 2.25f, 2.5f, 4.0f, 4.1f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(sve.i, 1, 10) {
        data, hidat0(1, sve.i);
      }
    }
    {
      static const float values[] = {
        2.0f, 3.0f, 5.0f, 6.0f, 7.0f, 8.0f, 8.0f, 10.0f, 10.0f, 10.0f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(sve.i, 1, 10) {
        data, hidat0(2, sve.i);
      }
    }
    {
      static const float values[] = {
        1.0f, 0.8f, 0.8f, 0.7f, 0.45f, 0.215f, 0.21f, 0.19f, 0.19f, 0.19f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(sve.i, 1, 10) {
        data, hidat0(3, sve.i);
      }
    }
    {
      static const float values[] = {
        0.35f, 0.35f, 0.3f, 0.3f, 0.3f, 0.3f, 0.5f, 0.6f, 0.6f, 0.6f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(sve.i, 1, 10) {
        data, hidat0(4, sve.i);
      }
    }
    {
      static const float values[] = {
        23.8f, 24.0f, 26.0f, 26.2f, 27.0f, 28.5f, 28.5f, 28.5f, 28.5f, 28.5f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(sve.i, 1, 10) {
        data, hidat0(5, sve.i);
      }
    }
    {
      fem::data_values data((values, 40*datum(0.0f)));
      FEM_DO_SAFE(sve.j, 6, 9) {
        FEM_DO_SAFE(sve.i, 1, 10) {
          data, hidat0(sve.j, sve.i);
        }
      }
    }
    {
      static const float values[] = {
        5.0f, 20.0f, 53.0f, 62.0f, 100.0f, 200.0f, 546.0f, 900.0f,
          1800.0f, 4000.0f
      };
      fem::data_of_type<float> data(FEM_VALUES_AND_SIZE);
      FEM_DO_SAFE(sve.i, 1, 10) {
        data, hidat0(10, sve.i);
      }
    }
    fem::data((values, 10*datum(0.0f))), hidat;
  }
  //Cc      SAVE /BVEG1/
  //Cc      SAVE /SEDVAX/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /HMAIN1/
  //Cc      SAVE /HMAIN2/
  //Cc      SAVE /HSTRNG/
  //Cc      SAVE /hjcrdn/
  //Cc      SAVE /HJJET1/
  //Cc      SAVE /HJJET2/
  //Cc      SAVE /HIJDAT/
  //Cc      SAVE /HPINT/
  //C...give all the switchs and parameters the default values
  //Clin-4/2008 input.ampt provides NSEED for AMPT:
  //C        DATA NSEED/74769375/
  //C
  //C...initialize all the data common blocks
  //Clin-4/26/01
  //C        DATA KATT/520000*0/PATT/520000*0.0/
  //C
  //Clin-4/2008
  //C        DATA NSG/0/,NJSG/900*0/,IASG/2700*0/,K1SG/90000*0/,K2SG/90000*0/
  //C     &       ,PXSG/90000*0.0/,PYSG/90000*0.0/,PZSG/90000*0.0/
  //C     &       ,PESG/90000*0.0/,PMSG/90000*0.0/
}

} // namespace AMPT
