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
iarflv(...)
{
  throw std::runtime_error(
    "Missing function implementation: iarflv");
}

void
idint(...)
{
  throw std::runtime_error(
    "Missing function implementation: idint");
}

void
index1(...)
{
  throw std::runtime_error(
    "Missing function implementation: index1");
}

void
invflv(...)
{
  throw std::runtime_error(
    "Missing function implementation: invflv");
}

void
lorenz(...)
{
  throw std::runtime_error(
    "Missing function implementation: lorenz");
}

void
ludecy(...)
{
  throw std::runtime_error(
    "Missing function implementation: ludecy");
}

void
lulist(...)
{
  throw std::runtime_error(
    "Missing function implementation: lulist");
}

void
ranart(...)
{
  throw std::runtime_error(
    "Missing function implementation: ranart");
}

void
rotate(...)
{
  throw std::runtime_error(
    "Missing function implementation: rotate");
}

void
ulmass(...)
{
  throw std::runtime_error(
    "Missing function implementation: ulmass");
}

struct common_snn
{
  float efrm;
  int npart1;
  int npart2;
  float epsipz;
  float epsipt;
  float pzproj;
  float pztarg;

  common_snn() :
    efrm(fem::float0),
    npart1(fem::int0),
    npart2(fem::int0),
    epsipz(fem::float0),
    epsipt(fem::float0),
    pzproj(fem::float0),
    pztarg(fem::float0)
  {}
};

struct common_lastt
{
  int itimeh;
  float bimp;

  common_lastt() :
    itimeh(fem::int0),
    bimp(fem::float0)
  {}
};

struct common_hbt
{
  static const int maxstr = 150001;

  arr<int> lblast;
  arr<float, 2> xlast;
  arr<float, 2> plast;
  int nlast;

  common_hbt() :
    lblast(dimension(maxstr), fem::fill0),
    xlast(dimension(4, maxstr), fem::fill0),
    plast(dimension(4, maxstr), fem::fill0),
    nlast(fem::int0)
  {}
};

const int common_hbt::maxstr;

struct common_oscar1
{
  int iap;
  int izp;
  int iat;
  int izt;

  common_oscar1() :
    iap(fem::int0),
    izp(fem::int0),
    iat(fem::int0),
    izt(fem::int0)
  {}
};

struct common_oscar2
{
  fem::str<8> frame;
  fem::str<25> amptvn;

  common_oscar2() :
    frame(fem::char0),
    amptvn(fem::char0)
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

struct common_input1
{
  int masspr;
  int massta;
  int iseed;
  int iavoid;
  float dt;

  common_input1() :
    masspr(fem::int0),
    massta(fem::int0),
    iseed(fem::int0),
    iavoid(fem::int0),
    dt(fem::float0)
  {}
};

struct common_aa
{
  static const int maxstr = 150001;

  arr<float, 2> r;

  common_aa() :
    r(dimension(3, maxstr), fem::fill0)
  {}
};

const int common_aa::maxstr;

struct common_bb
{
  static const int maxstr = 150001;

  arr<float, 2> p;

  common_bb() :
    p(dimension(3, maxstr), fem::fill0)
  {}
};

const int common_bb::maxstr;

struct common_cc
{
  static const int maxstr = 150001;

  arr<float> e;

  common_cc() :
    e(dimension(maxstr), fem::fill0)
  {}
};

const int common_cc::maxstr;

struct common_ee
{
  static const int maxstr = 150001;

  arr<int> id;
  arr<int> lb;

  common_ee() :
    id(dimension(maxstr), fem::fill0),
    lb(dimension(maxstr), fem::fill0)
  {}
};

const int common_ee::maxstr;

struct common_tdecay
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float> tfdcy;
  arr<float, 2> tfdpi;
  arr<float> tft;

  common_tdecay() :
    tfdcy(dimension(maxstr), fem::fill0),
    tfdpi(dimension(maxstr, maxr), fem::fill0),
    tft(dimension(maxstr), fem::fill0)
  {}
};

const int common_tdecay::maxstr;
const int common_tdecay::maxr;

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

struct common_ftmax
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float> ftsv;
  arr<float, 2> ftsvt;

  common_ftmax() :
    ftsv(dimension(maxstr), fem::fill0),
    ftsvt(dimension(maxstr, maxr), fem::fill0)
  {}
};

const int common_ftmax::maxstr;
const int common_ftmax::maxr;

struct common_dpert
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 2> dpertt;
  arr<float> dpertp;
  arr<float> dplast;
  arr<float> dpdcy;
  arr<float, 2> dpdpi;
  arr<float, 2> dpt;
  arr<float, 2> dpp1;
  arr<float, 2> dppion;

  common_dpert() :
    dpertt(dimension(maxstr, maxr), fem::fill0),
    dpertp(dimension(maxstr), fem::fill0),
    dplast(dimension(maxstr), fem::fill0),
    dpdcy(dimension(maxstr), fem::fill0),
    dpdpi(dimension(maxstr, maxr), fem::fill0),
    dpt(dimension(maxstr, maxr), fem::fill0),
    dpp1(dimension(maxstr, maxr), fem::fill0),
    dppion(dimension(maxstr, maxr), fem::fill0)
  {}
};

const int common_dpert::maxstr;
const int common_dpert::maxr;

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

struct common_para8
{
  int idpert;
  int npertd;
  int idxsec;

  common_para8() :
    idpert(fem::int0),
    npertd(fem::int0),
    idxsec(fem::int0)
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

struct common_lor
{
  double enenew;
  double pxnew;
  double pynew;
  double pznew;

  common_lor() :
    enenew(fem::double0),
    pxnew(fem::double0),
    pynew(fem::double0),
    pznew(fem::double0)
  {}
};

struct common_decom
{
  arr<double, 2> ptwo;

  common_decom() :
    ptwo(dimension(2, 5), fem::fill0)
  {}
};

struct common_rndf77
{
  int nseed;

  common_rndf77() :
    nseed(fem::int0)
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

struct common_embed
{
  int iembed;
  int nsembd;
  float pxqembd;
  float pyqembd;
  float xembd;
  float yembd;
  float psembd;
  float tmaxembd;
  float phidecomp;

  common_embed() :
    iembed(fem::int0),
    nsembd(fem::int0),
    pxqembd(fem::float0),
    pyqembd(fem::float0),
    xembd(fem::float0),
    yembd(fem::float0),
    psembd(fem::float0),
    tmaxembd(fem::float0),
    phidecomp(fem::float0)
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

struct common_loclco
{
  arr<double> gxp;
  arr<double> gyp;
  arr<double> gzp;
  arr<double> ftp;
  arr<double> pxp;
  arr<double> pyp;
  arr<double> pzp;
  arr<double> pep;
  arr<double> pmp;

  common_loclco() :
    gxp(dimension(3), fem::fill0),
    gyp(dimension(3), fem::fill0),
    gzp(dimension(3), fem::fill0),
    ftp(dimension(3), fem::fill0),
    pxp(dimension(3), fem::fill0),
    pyp(dimension(3), fem::fill0),
    pzp(dimension(3), fem::fill0),
    pep(dimension(3), fem::fill0),
    pmp(dimension(3), fem::fill0)
  {}
};

struct common_prtn23
{
  arr<double> gxp0;
  arr<double> gyp0;
  arr<double> gzp0;
  double ft0fom;

  common_prtn23() :
    gxp0(dimension(3), fem::fill0),
    gyp0(dimension(3), fem::fill0),
    gzp0(dimension(3), fem::fill0),
    ft0fom(fem::double0)
  {}
};

struct common_coal
{
  double dpcoal;
  double drcoal;
  double ecritl;

  common_coal() :
    dpcoal(fem::double0),
    drcoal(fem::double0),
    ecritl(fem::double0)
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

struct common_arprnt
{
  arr<float> arpar1;
  arr<int> iapar2;
  arr<float> arint1;
  arr<int> iaint2;

  common_arprnt() :
    arpar1(dimension(100), fem::fill0),
    iapar2(dimension(50), fem::fill0),
    arint1(dimension(100), fem::fill0),
    iaint2(dimension(50), fem::fill0)
  {}
};

struct common_nzpc
{
  int nattzp;

  common_nzpc() :
    nattzp(fem::int0)
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

struct common_input2
{
  int ilab;
  int manyb;
  int ntmax;
  int icoll;
  int insys;
  int ipot;
  int mode;
  int imomen;
  int nfreq;
  int icflow;
  int icrho;
  int icou;
  int kpoten;
  int kmul;

  common_input2() :
    ilab(fem::int0),
    manyb(fem::int0),
    ntmax(fem::int0),
    icoll(fem::int0),
    insys(fem::int0),
    ipot(fem::int0),
    mode(fem::int0),
    imomen(fem::int0),
    nfreq(fem::int0),
    icflow(fem::int0),
    icrho(fem::int0),
    icou(fem::int0),
    kpoten(fem::int0),
    kmul(fem::int0)
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

struct common_pa
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 3> rpion;

  common_pa() :
    rpion(dimension(3, maxstr, maxr), fem::fill0)
  {}
};

const int common_pa::maxstr;
const int common_pa::maxr;

struct common_pb
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 3> ppion;

  common_pb() :
    ppion(dimension(3, maxstr, maxr), fem::fill0)
  {}
};

const int common_pb::maxstr;
const int common_pb::maxr;

struct common_pc
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 2> epion;

  common_pc() :
    epion(dimension(maxstr, maxr), fem::fill0)
  {}
};

const int common_pc::maxstr;
const int common_pc::maxr;

struct common_pd
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<int, 2> lpion;

  common_pd() :
    lpion(dimension(maxstr, maxr), fem::fill0)
  {}
};

const int common_pd::maxstr;
const int common_pd::maxr;

struct common_resdcy
{
  int nsav;
  int iksdcy;

  common_resdcy() :
    nsav(fem::int0),
    iksdcy(fem::int0)
  {}
};

struct common_leadng
{
  int lb1;
  float px1;
  float py1;
  float pz1;
  float em1;
  float e1;
  float xfnl;
  float yfnl;
  float zfnl;
  float tfnl;
  float px1n;
  float py1n;
  float pz1n;
  float dp1n;

  common_leadng() :
    lb1(fem::int0),
    px1(fem::float0),
    py1(fem::float0),
    pz1(fem::float0),
    em1(fem::float0),
    e1(fem::float0),
    xfnl(fem::float0),
    yfnl(fem::float0),
    zfnl(fem::float0),
    tfnl(fem::float0),
    px1n(fem::float0),
    py1n(fem::float0),
    pz1n(fem::float0),
    dp1n(fem::float0)
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

struct common_frzprc
{
  static const int maxptn = 400001;

  arr<double> gxfrz;
  arr<double> gyfrz;
  arr<double> gzfrz;
  arr<double> ftfrz;
  arr<double> pxfrz;
  arr<double> pyfrz;
  arr<double> pzfrz;
  arr<double> efrz;
  arr<double> xmfrz;
  arr<double> tfrz;
  arr<int> ifrz;
  arr<int> idfrz;
  int itlast;

  common_frzprc() :
    gxfrz(dimension(maxptn), fem::fill0),
    gyfrz(dimension(maxptn), fem::fill0),
    gzfrz(dimension(maxptn), fem::fill0),
    ftfrz(dimension(maxptn), fem::fill0),
    pxfrz(dimension(maxptn), fem::fill0),
    pyfrz(dimension(maxptn), fem::fill0),
    pzfrz(dimension(maxptn), fem::fill0),
    efrz(dimension(maxptn), fem::fill0),
    xmfrz(dimension(maxptn), fem::fill0),
    tfrz(dimension(302), fem::fill0),
    ifrz(dimension(maxptn), fem::fill0),
    idfrz(dimension(maxptn), fem::fill0),
    itlast(fem::int0)
  {}
};

const int common_frzprc::maxptn;

struct common_prec4
{
  static const int maxptn = 400001;

  arr<double> vx;
  arr<double> vy;
  arr<double> vz;

  common_prec4() :
    vx(dimension(maxptn), fem::fill0),
    vy(dimension(maxptn), fem::fill0),
    vz(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec4::maxptn;

struct common_prec5
{
  static const int maxptn = 400001;

  arr<double> eta;
  arr<double> rap;
  arr<double> tau;

  common_prec5() :
    eta(dimension(maxptn), fem::fill0),
    rap(dimension(maxptn), fem::fill0),
    tau(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec5::maxptn;

struct common_ilist5
{
  static const int maxptn = 400001;

  arr<double> ct;
  arr<double> ot;
  double tlarge;

  common_ilist5() :
    ct(dimension(maxptn), fem::fill0),
    ot(dimension(maxptn), fem::fill0),
    tlarge(fem::double0)
  {}
};

const int common_ilist5::maxptn;

struct common_hflow
{
  arr<double, 2> v2h;
  arr<double, 2> xnhadr;
  arr<double, 2> eth;
  arr<double, 2> v2h2;
  arr<double, 2> s2h;

  common_hflow() :
    v2h(dimension(30, 3), fem::fill0),
    xnhadr(dimension(30, 3), fem::fill0),
    eth(dimension(30, 3), fem::fill0),
    v2h2(dimension(30, 3), fem::fill0),
    s2h(dimension(30, 3), fem::fill0)
  {}
};

struct common_ebe
{
  arr<double> v2hp;
  arr<double> xnhadp;
  arr<double> v2hsum;
  arr<double> v2h2sm;

  common_ebe() :
    v2hp(dimension(3), fem::fill0),
    xnhadp(dimension(3), fem::fill0),
    v2hsum(dimension(3), fem::fill0),
    v2h2sm(dimension(3), fem::fill0)
  {}
};

struct common_run
{
  int num;

  common_run() :
    num(fem::int0)
  {}
};

struct common_rr
{
  static const int maxr = 1;

  arr<int> massr;

  common_rr() :
    massr(dim1(0, maxr), fem::fill0)
  {}
};

const int common_rr::maxr;

struct common_arerc1
{
  static const int maxr = 1;

  arr<int> multi1;

  common_arerc1() :
    multi1(dimension(maxr), fem::fill0)
  {}
};

const int common_arerc1::maxr;

struct common_arprc1
{
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<int, 2> ityp1;
  arr<float, 2> gx1;
  arr<float, 2> gy1;
  arr<float, 2> gz1;
  arr<float, 2> ft1;
  arr<float, 2> px1;
  arr<float, 2> py1;
  arr<float, 2> pz1;
  arr<float, 2> ee1;
  arr<float, 2> xm1;

  common_arprc1() :
    ityp1(dimension(maxstr, maxr), fem::fill0),
    gx1(dimension(maxstr, maxr), fem::fill0),
    gy1(dimension(maxstr, maxr), fem::fill0),
    gz1(dimension(maxstr, maxr), fem::fill0),
    ft1(dimension(maxstr, maxr), fem::fill0),
    px1(dimension(maxstr, maxr), fem::fill0),
    py1(dimension(maxstr, maxr), fem::fill0),
    pz1(dimension(maxstr, maxr), fem::fill0),
    ee1(dimension(maxstr, maxr), fem::fill0),
    xm1(dimension(maxstr, maxr), fem::fill0)
  {}
};

const int common_arprc1::maxstr;
const int common_arprc1::maxr;

struct common_iflow
{
  double v2i;
  double eti;
  double xmulti;
  double v2mi;
  double s2mi;
  double xmmult;
  double v2bi;
  double s2bi;
  double xbmult;

  common_iflow() :
    v2i(fem::double0),
    eti(fem::double0),
    xmulti(fem::double0),
    v2mi(fem::double0),
    s2mi(fem::double0),
    xmmult(fem::double0),
    v2bi(fem::double0),
    s2bi(fem::double0),
    xbmult(fem::double0)
  {}
};

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

struct common_xyembed
{
  static const int nxymax = 10001;

  int nxyjet;
  arr<float, 2> xyjet;

  common_xyembed() :
    nxyjet(fem::int0),
    xyjet(dimension(nxymax, 2), fem::fill0)
  {}
};

const int common_xyembed::nxymax;

struct common :
  fem::common,
  common_snn,
  common_lastt,
  common_hbt,
  common_oscar1,
  common_oscar2,
  common_para7,
  common_input1,
  common_aa,
  common_bb,
  common_cc,
  common_ee,
  common_tdecay,
  common_arevt,
  common_hjglbr,
  common_ftmax,
  common_dpert,
  common_hparnt,
  common_para8,
  common_phihj,
  common_lor,
  common_decom,
  common_rndf77,
  common_hmain1,
  common_embed,
  common_hmain2,
  common_para1,
  common_prec1,
  common_ilist7,
  common_arprc,
  common_noprec,
  common_soft,
  common_anim,
  common_precpa,
  common_loclco,
  common_prtn23,
  common_coal,
  common_hjjet2,
  common_arprnt,
  common_nzpc,
  common_ludat1,
  common_input2,
  common_lujets,
  common_ludat2,
  common_ludat3,
  common_pa,
  common_pb,
  common_pc,
  common_pd,
  common_resdcy,
  common_leadng,
  common_phidcy,
  common_prec2,
  common_frzprc,
  common_prec4,
  common_prec5,
  common_ilist5,
  common_hflow,
  common_ebe,
  common_run,
  common_rr,
  common_arerc1,
  common_arprc1,
  common_iflow,
  common_frzout,
  common_hjcrdn,
  common_hjjet1,
  common_xyembed
{
  fem::cmn_sve hoscar_sve;
  fem::cmn_sve hbtout_sve;
  fem::cmn_sve decomp_sve;
  fem::cmn_sve htop_sve;
  fem::cmn_sve resmass_sve;
  fem::cmn_sve exchge_sve;
  fem::cmn_sve locldr_sve;
  fem::cmn_sve coales_sve;
  fem::cmn_sve ptoh_sve;
  fem::cmn_sve getnp_sve;
  fem::cmn_sve resdec_sve;
  fem::cmn_sve local_sve;
  fem::cmn_sve inifrz_sve;
  fem::cmn_sve flowh_sve;
  fem::cmn_sve flowh0_sve;
  fem::cmn_sve iniflw_sve;
  fem::cmn_sve frztm_sve;
  fem::cmn_sve minijet_out_sve;
  fem::cmn_sve embedhighpt_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct hoscar_save
{
  fem::str<8> code;
  float ebeam;
  float ene;
  int i;
  int ievent;
  int nff;
  int ntestp;
  float phi;
  fem::str<8> reffra;
  float xmp;
  float xmt;

  hoscar_save() :
    code(fem::char0),
    ebeam(fem::float0),
    ene(fem::float0),
    i(fem::int0),
    ievent(fem::int0),
    nff(fem::int0),
    ntestp(fem::int0),
    phi(fem::float0),
    reffra(fem::char0),
    xmp(fem::float0),
    xmt(fem::float0)
  {}
};

//C
//C=======================================================================
void
hoscar(
  common& cmn)
{
  FEM_CMN_SVE(hoscar);
  common_write write(cmn);
  // COMMON snn
  float& efrm = cmn.efrm;
  // COMMON hbt
  const int maxstr = 150001;
  arr_cref<int> lblast(cmn.lblast, dimension(maxstr));
  arr_cref<float, 2> xlast(cmn.xlast, dimension(4, maxstr));
  arr_cref<float, 2> plast(cmn.plast, dimension(4, maxstr));
  int& nlast = cmn.nlast;
  // COMMON oscar1
  int& iap = cmn.iap;
  int& izp = cmn.izp;
  int& iat = cmn.iat;
  int& izt = cmn.izt;
  // COMMON oscar2
  fem::str<8>& frame = cmn.frame;
  //
  // SAVE
  fem::str<8>& code = sve.code;
  float& ebeam = sve.ebeam;
  float& ene = sve.ene;
  int& i = sve.i;
  int& ievent = sve.ievent;
  int& nff = sve.nff;
  int& ntestp = sve.ntestp;
  float& phi = sve.phi;
  fem::str<8>& reffra = sve.reffra;
  float& xmp = sve.xmp;
  float& xmt = sve.xmt;
  //
  if (is_called_first_time) {
    nff = 0;
  }
  //C
  //Cc      SAVE /snn/
  //Cc      SAVE /lastt/
  //Cc      SAVE /hbt/
  //Cc      SAVE /oscar1/
  //Cc      SAVE /oscar2/
  //C
  //C       file header
  const float amp = 0.93828f;
  const float amn = 0.939457f;
  if (nff == 0) {
    write(19, "(a10)"), "OSCAR1997A";
    write(19, "(a12)"), "final_id_p_x";
    code = "AMPT";
    if (frame == "CMS") {
      reffra = "nncm";
      xmp = (amp * izp + amn * (iap - izp)) / iap;
      xmt = (amp * izt + amn * (iat - izt)) / iat;
      ebeam = (fem::pow2(efrm) - fem::pow2(xmp) - fem::pow2(xmt)) / 2.f / xmt;
    }
    else if (frame == "LAB") {
      reffra = "lab";
      ebeam = efrm;
    }
    else {
      reffra = "unknown";
      ebeam = 0.f;
    }
    ntestp = 1;
    write(19,
      "(a4,1x,a20,1x,'(',i3,',',i3,')+(',i3,',',i3,')',2x,a4,2x,e10.4,2x,i8)"),
      code, cmn.amptvn, iap, izp, iat, izt, reffra, ebeam, ntestp;
    nff = 1;
    ievent = 1;
    phi = 0.f;
    if (frame == "CMS") {
      write(19,
        "('# Center-of-mass energy/nucleon-pair is',f12.3,'GeV')"),
        efrm;
    }
  }
  //C       comment
  //C       event header
  write(19, "(i10,2x,i10,2x,f8.3,2x,f8.3)"), ievent, nlast, cmn.bimp, phi;
  //C       particles
  FEM_DO_SAFE(i, 1, nlast) {
    ene = fem::sqrt(fem::pow2(plast(1, i)) + fem::pow2(plast(2, i)) +
      fem::pow2(plast(3, i)) + fem::pow2(plast(4, i)));
    write(19, "(i10,2x,i10,2x,9(e12.6,2x))"), i, invflv(lblast(i)),
      plast(1, i), plast(2, i), plast(3, i), ene, plast(4, i), xlast(1,
      i), xlast(2, i), xlast(3, i), xlast(4, i);
  }
  ievent++;
  //C
}

struct hbtout_save
{
  static const int maxstr = 150001;

  float deltat;
  float dr;
  float ene;
  int i;
  int ii;
  int ip;
  int ip2;
  int iplast;
  arr<int> lastkp;
  int ndpert;
  arr<int> newkp;
  arr<float> xnew;

  hbtout_save() :
    deltat(fem::float0),
    dr(fem::float0),
    ene(fem::float0),
    i(fem::int0),
    ii(fem::int0),
    ip(fem::int0),
    ip2(fem::int0),
    iplast(fem::int0),
    lastkp(dimension(maxstr), fem::fill0),
    ndpert(fem::int0),
    newkp(dimension(maxstr), fem::fill0),
    xnew(dimension(3), fem::fill0)
  {}
};

const int hbtout_save::maxstr;

//C.................... linana.f
//C=======================================================================
//C     10/26/01 update freezeout positions in case of interactions:
//Clin-3/2009 Note: freezeout spacetime values cannot be trusted for K0S & K0L
//C     as K0S/K0L are converted from K+/K- by hand at the end of hadron cascade.
void
hbtout(
  common& cmn,
  int const& nnew,
  int const& nt,
  int const& ntmax)
{
  FEM_CMN_SVE(hbtout);
  common_write write(cmn);
  const int maxstr = 150001;
  arr_ref<int> lblast(cmn.lblast, dimension(maxstr));
  arr_ref<float, 2> xlast(cmn.xlast, dimension(4, maxstr));
  arr_ref<float, 2> plast(cmn.plast, dimension(4, maxstr));
  int& nlast = cmn.nlast;
  float& dt = cmn.dt;
  arr_cref<float, 2> r(cmn.r, dimension(3, maxstr));
  arr_cref<float, 2> p(static_cast<common_bb&>(cmn).p, dimension(3, maxstr));
  arr_cref<float> e(cmn.e, dimension(maxstr));
  arr_cref<int> lb(cmn.lb, dimension(maxstr));
  float& bimp = cmn.bimp;
  arr_cref<float> tfdcy(cmn.tfdcy, dimension(maxstr));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  int& npart1 = cmn.npart1;
  int& npart2 = cmn.npart2;
  float& epsipz = cmn.epsipz;
  float& epsipt = cmn.epsipt;
  int& nelt = cmn.nelt;
  int& ninthj = cmn.ninthj;
  int& nelp = cmn.nelp;
  int& ninp = cmn.ninp;
  arr_cref<float> ftsv(cmn.ftsv, dimension(maxstr));
  arr_cref<float> dpertp(cmn.dpertp, dimension(maxstr));
  arr_ref<float> dplast(cmn.dplast, dimension(maxstr));
  int& idpert = cmn.idpert;
  //
  float& deltat = sve.deltat;
  float& dr = sve.dr;
  float& ene = sve.ene;
  int& i = sve.i;
  int& ii = sve.ii;
  int& ip = sve.ip;
  int& ip2 = sve.ip2;
  int& iplast = sve.iplast;
  arr_ref<int> lastkp(sve.lastkp, dimension(maxstr));
  int& ndpert = sve.ndpert;
  arr_ref<int> newkp(sve.newkp, dimension(maxstr));
  arr_ref<float> xnew(sve.xnew, dimension(3));
  const float oneminus = 0.99999f;
  const float oneplus = 1.00001f;
  static const char* format_200 = "(i6,2(1x,f8.3),1x,f11.4,1x,f6.3,4(1x,f8.2))";
  static const char* format_250 =
    "(i5,2(1x,f8.3),1x,f10.3,2(1x,f7.1),1x,f8.2,1x,f7.2,1x,e10.4)";
  //C
  //Clin-5/2008 give tolerance to regular particles (perturbative probability 1):
  //Cc      SAVE /para7/
  //Cc      SAVE /hbt/
  //Cc      SAVE /input1/
  //Cc      SAVE /AA/
  //Cc      SAVE /BB/
  //Cc      SAVE /CC/
  //Cc      SAVE /EE/
  //Cc      SAVE /lastt/
  //Cc      SAVE /tdecay/
  //Cc      SAVE /AREVT/
  //Cc      SAVE /snn/
  //Cc      SAVE /HJGLBR/
  //Clin-12/14/03:
  //Clin-2/2012:
  //C
  FEM_DO_SAFE(i, 1, fem::max0(nlast, nnew)) {
    lastkp(i) = 0;
  }
  FEM_DO_SAFE(i, 1, nnew) {
    newkp(i) = 0;
  }
  //C     for each of the particles, search the freezeout record (common /hbt/)
  //C     to find & keep those which do not have interactions during this timestep:
  FEM_DO_SAFE(ip, 1, nnew) {
    FEM_DO_SAFE(iplast, 1, nlast) {
      if (p(1, ip) == plast(1, iplast) && p(2, ip) == plast(2,
          iplast) && p(3, ip) == plast(3, iplast) && e(ip) == plast(4,
          iplast) && lb(ip) == lblast(iplast) && dpertp(ip) == dplast(
          iplast) && lastkp(iplast) == 0) {
        //Clin-5/2008 modified below to the above in case we have perturbative particles:
        //C     5           lastkp(iplast).eq.0) then
        deltat = nt * dt - xlast(4, iplast);
        ene = fem::sqrt(fem::pow2(plast(1, iplast)) + fem::pow2(plast(2,
          iplast)) + fem::pow2(plast(3, iplast)) + fem::pow2(plast(4,
          iplast)));
        //C     xnew gives the coordinate if a particle free-streams to current time:
        FEM_DO_SAFE(ii, 1, 3) {
          xnew(ii) = xlast(ii, iplast) + plast(ii, iplast) / ene * deltat;
        }
        dr = fem::sqrt(fem::pow2((r(1, ip) - xnew(1))) + fem::pow2((r(2,
          ip) - xnew(2))) + fem::pow2((r(3, ip) - xnew(3))));
        //C     find particles with dp=0 and dr<0.01, considered to be those
        //C     without any interactions during this timestep,
        //C     thus keep their last positions and time:
        if (dr <= 0.01f) {
          lastkp(iplast) = 1;
          newkp(ip) = 1;
          //C                  if(lb(ip).eq.41) then
          //C                write(95,*) 'nt,ip,px,x=',nt,ip,p(1,ip),r(1,ip),ftsv(ip)
          //C                write(95,*) 'xnew=',xnew(1),xnew(2),xnew(3),xlast(4,ip)
          //C                  endif
          //Clin-5/2009 Take care of formation time of particles read in at nt=ntmax-1:
          if (nt == ntmax && ftsv(ip) > ((ntmax - 1) * dt)) {
            xlast(4, iplast) = ftsv(ip);
          }
          goto statement_100;
        }
      }
    }
    statement_100:;
  }
  //C     for current particles with interactions, fill their current info in
  //C     the freezeout record (if that record entry needs not to be kept):
  FEM_DO_SAFE(ip, 1, nnew) {
    if (newkp(ip) == 0) {
      FEM_DO_SAFE(iplast, 1, nnew) {
        if (lastkp(iplast) == 0) {
          //Ctest off: write collision info
          //C                  if(lb(ip).eq.41) then
          //C                     write(95,*) 'nt,lb(ip)=',nt,lb(ip)
          //C                  write(95,*) '  last p=',plast(1,iplast),
          //C     1 plast(2,iplast),plast(3,iplast),plast(4,iplast)
          //C                  write(95,*) '  after p=',p(1,ip),p(2,ip),p(3,ip),e(ip)
          //C                  write(95,*) 'after x=',r(1,ip),r(2,ip),r(3,ip),ftsv(ip)
          //C                  endif
          //C
          xlast(1, iplast) = r(1, ip);
          xlast(2, iplast) = r(2, ip);
          xlast(3, iplast) = r(3, ip);
          xlast(4, iplast) = nt * dt;
          //C
          if (nt == ntmax) {
            //C     freezeout time for decay daughters at the last timestep
            //C     needs to include the decay time of the parent:
            if (tfdcy(ip) > (ntmax * dt + 0.001f)) {
              xlast(4, iplast) = tfdcy(ip);
              //C     freezeout time for particles unformed at the next-to-last timestep
              //C     needs to be their formation time instead of (ntmax*dt):
            }
            else if (ftsv(ip) > ((ntmax - 1) * dt)) {
              xlast(4, iplast) = ftsv(ip);
            }
          }
          plast(1, iplast) = p(1, ip);
          plast(2, iplast) = p(2, ip);
          plast(3, iplast) = p(3, ip);
          plast(4, iplast) = e(ip);
          lblast(iplast) = lb(ip);
          lastkp(iplast) = 1;
          //Clin-5/2008:
          dplast(iplast) = dpertp(ip);
          goto statement_150;
        }
      }
    }
    statement_150:;
  }
  //C     if the current particle list is shorter than the freezeout record,
  //C     condense the last-collision record by filling new record from 1 to nnew,
  //C     and label these entries as keep:
  if (nnew < nlast) {
    FEM_DO_SAFE(iplast, 1, nlast) {
      if (lastkp(iplast) == 0) {
        FEM_DO_SAFE(ip2, iplast + 1, nlast) {
          if (lastkp(ip2) == 1) {
            xlast(1, iplast) = xlast(1, ip2);
            xlast(2, iplast) = xlast(2, ip2);
            xlast(3, iplast) = xlast(3, ip2);
            xlast(4, iplast) = xlast(4, ip2);
            plast(1, iplast) = plast(1, ip2);
            plast(2, iplast) = plast(2, ip2);
            plast(3, iplast) = plast(3, ip2);
            plast(4, iplast) = plast(4, ip2);
            lblast(iplast) = lblast(ip2);
            lastkp(iplast) = 1;
            //Clin-5/2008:
            dplast(iplast) = dplast(ip2);
            goto statement_170;
          }
        }
      }
      statement_170:;
    }
  }
  nlast = nnew;
  //Ctest off look inside each NT timestep (for debugging purpose):
  //C      do ip=1,nlast
  //C         write(99,*) ' p ',nt,ip,lblast(ip),plast(1,ip),
  //C     1        plast(2,ip),plast(3,ip),plast(4,ip),dplast(ip)
  //C         write(99,*) '  x ',nt,ip,lblast(ip),xlast(1,ip),
  //C     1        xlast(2,ip),xlast(3,ip),xlast(4,ip),dplast(ip)
  //C      enddo
  //C
  if (nt == ntmax) {
    //Clin-5/2008 find final number of perturbative particles (deuterons only):
    ndpert = 0;
    FEM_DO_SAFE(ip, 1, nlast) {
      if (dplast(ip) > oneminus && dplast(ip) < oneplus) {
      }
      else {
        ndpert++;
      }
    }
    //C
    //C         write(16,190) IAEVT,IARUN,nlast,bimp,npart1,npart2,
    //C     1 NELP,NINP,NELT,NINTHJ
    //Clin-2/2012:
    //C         write(16,190) IAEVT,IARUN,nlast-ndpert,bimp,npart1,npart2,
    //C     1 NELP,NINP,NELT,NINTHJ
    write(16, "(3(i7),f10.4,5x,6(i4),5x,f7.4)"), iaevt, iarun,
      nlast - ndpert, bimp, npart1, npart2, nelp, ninp, nelt, ninthj,
      cmn.phirp;
    //Clin-5/2008 write out perturbatively-produced particles (deuterons only):
    if (idpert == 1 || idpert == 2) {
      write(90, "(3(i7),f10.4,5x,6(i4))"), iaevt, iarun, ndpert,
        bimp, npart1, npart2, nelp, ninp, nelt, ninthj;
    }
    FEM_DO_SAFE(ip, 1, nlast) {
      //Clin-12/14/03   No formation time for spectator projectile or target nucleons,
      //C     see ARINI1 in 'amptsub.f':
      //Clin-3/2009 To be consistent with new particles produced in hadron cascade
      //C     that are limited by the time-resolution (DT) of the hadron cascade,
      //C     freezeout time of spectator projectile or target nucleons is written as
      //C     DT as they are read at the 1st timestep and then propagated to time DT:
      //C
      //Clin-9/2011 determine spectator nucleons consistently
      //C            if(plast(1,ip).eq.0.and.plast(2,ip).eq.0
      //C     1           .and.(sqrt(plast(3,ip)**2+plast(4,ip)**2)*2/HINT1(1))
      //C     2           .gt.0.99.and.(lblast(ip).eq.1.or.lblast(ip).eq.2)) then
      if (fem::abs(plast(1, ip)) <= epsipt && fem::abs(plast(2,
          ip)) <= epsipt && (plast(3, ip) > fem::amax1(0.f,
          cmn.pzproj - epsipz) || plast(3, ip) < (-cmn.pztarg +
          epsipz)) && (lblast(ip) == 1 || lblast(ip) == 2)) {
        //Clin-5/2008 perturbatively-produced particles (currently only deuterons)
        //C     are written to ana/ampt_pert.dat (without the column for the mass);
        //C     ana/ampt.dat has regularly-produced particles (including deuterons);
        //C     these two sets of deuteron data are close to each other(but not the same
        //C     because of the bias from triggering the perturbative production);
        //C     ONLY use one data set for analysis to avoid double-counting:
        if (dplast(ip) > oneminus && dplast(ip) < oneplus) {
          write(16, format_200), invflv(lblast(ip)), plast(1, ip),
            plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip),
            xlast(2, ip), xlast(3, ip), xlast(4, ip);
          //Clin-12/14/03-end
        }
        else {
          if (idpert == 1 || idpert == 2) {
            write(90, format_250), invflv(lblast(ip)), plast(1, ip),
              plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip),
              xlast(3, ip), xlast(4, ip);
          }
          else {
            write(99, star), "Unexpected perturbative particles";
          }
        }
      }
      else if (fem::amax1(fem::abs(xlast(1, ip)), fem::abs(xlast(2,
        ip)), fem::abs(xlast(3, ip)), fem::abs(xlast(4,
        ip))) < 9999) {
        if (dplast(ip) > oneminus && dplast(ip) < oneplus) {
          write(16, format_200), invflv(lblast(ip)), plast(1, ip),
            plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip),
            xlast(2, ip), xlast(3, ip), xlast(4, ip);
        }
        else {
          if (idpert == 1 || idpert == 2) {
            write(90, format_250), invflv(lblast(ip)), plast(1, ip),
              plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip),
              xlast(3, ip), xlast(4, ip), dplast(ip);
          }
          else {
            write(99, star), "Unexpected perturbative particles";
          }
        }
      }
      else {
        //C     change format for large numbers:
        if (dplast(ip) > oneminus && dplast(ip) < oneplus) {
          write(16, "(i6,2(1x,f8.3),1x,f11.4,1x,f6.3,4(1x,e8.2))"),
            invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3,
            ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3,
            ip), xlast(4, ip);
        }
        else {
          if (idpert == 1 || idpert == 2) {
            write(90,
              "(i5,2(1x,f8.3),1x,f10.3,4(1x,e8.2),1x,e10.4)"), invflv(
              lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip),
              xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip),
              dplast(ip);
          }
          else {
            write(99, star), "Unexpected perturbative particles";
          }
        }
      }
    }
    if (cmn.ioscar == 1) {
      hoscar(cmn);
    }
  }
  //Clin-3/2009 improve the output accuracy of Pz
  //C
}

struct decomp_save
{
  double beta2;
  double dbex;
  double dbey;
  double dbez;
  double dcth;
  double de0;
  double de1;
  double de2;
  double dpcm;
  double dphi;
  double dpx;
  double dpy;
  double dpz;
  double ds;
  double gam;

  decomp_save() :
    beta2(fem::double0),
    dbex(fem::double0),
    dbey(fem::double0),
    dbez(fem::double0),
    dcth(fem::double0),
    de0(fem::double0),
    de1(fem::double0),
    de2(fem::double0),
    dpcm(fem::double0),
    dphi(fem::double0),
    dpx(fem::double0),
    dpy(fem::double0),
    dpz(fem::double0),
    ds(fem::double0),
    gam(fem::double0)
  {}
};

//C
//C=======================================================================
void
decomp(
  common& cmn,
  double const& px0,
  double const& py0,
  double const& pz0,
  double const& xm0,
  int const& i,
  int const& /* itq1 */)
{
  FEM_CMN_SVE(decomp);
  common_write write(cmn);
  // COMMON lor
  double& enenew = cmn.enenew;
  double& pxnew = cmn.pxnew;
  double& pynew = cmn.pynew;
  double& pznew = cmn.pznew;
  // COMMON hparnt
  arr_cref<float> hipr1(cmn.hipr1, dimension(100));
  // COMMON decom
  arr_ref<double, 2> ptwo(cmn.ptwo, dimension(2, 5));
  // COMMON rndf77
  int& nseed = cmn.nseed;
  // COMMON hmain1
  int& natt = cmn.natt;
  // COMMON embed
  int& iembed = cmn.iembed;
  int& nsembd = cmn.nsembd;
  //
  // SAVE
  double& beta2 = sve.beta2;
  double& dbex = sve.dbex;
  double& dbey = sve.dbey;
  double& dbez = sve.dbez;
  double& dcth = sve.dcth;
  double& de0 = sve.de0;
  double& de1 = sve.de1;
  double& de2 = sve.de2;
  double& dpcm = sve.dpcm;
  double& dphi = sve.dphi;
  double& dpx = sve.dpx;
  double& dpy = sve.dpy;
  double& dpz = sve.dpz;
  double& ds = sve.ds;
  double& gam = sve.gam;
  //
  //C
  //Clin-8/2015 changed ptwo(2,5) and related variables to double precision
  //C     to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID or IEEE_OVERFLOW_FLAG:
  //Cc      SAVE /lor/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /decom/
  //Cc      SAVE /RNDF77/
  //C
  dcth = fem::dble(ranart(nseed)) * 2.e0 - 1.e0;
  dphi = fem::dble(ranart(nseed) * hipr1(40)) * 2.e0;
  //Clin-6/2009 Added if embedding a high-Pt quark pair after string melting:
  if (iembed >= 1 && iembed <= 4) {
    //C     Decompose the parent high-Pt pion to q and qbar with an internal momentum
    //C     parallel to the pion direction so that one parton has ~the same hight Pt
    //C     and the other parton has a very soft Pt:
    //C     Note: htop() decomposes a meson to q as it(1) followed by qbar as it(2):
    if (i == (natt - 2 * nsembd) || i == (natt - 2 * nsembd - 1)) {
      dcth = 0.e0;
      dphi = fem::dble(cmn.phidecomp);
    }
  }
  //C
  ds = fem::pow2(xm0);
  dpcm = fem::dsqrt((ds - fem::pow2((ptwo(1, 5) + ptwo(2, 5)))) * (
    ds - fem::pow2((ptwo(1, 5) - ptwo(2, 5)))) / ds / 4e0);
  dpz = dpcm * dcth;
  dpx = dpcm * fem::dsqrt(1.e0 - fem::pow2(dcth)) * fem::dcos(dphi);
  dpy = dpcm * fem::dsqrt(1.e0 - fem::pow2(dcth)) * fem::dsin(dphi);
  de1 = fem::dsqrt(fem::pow2(ptwo(1, 5)) + fem::pow2(dpcm));
  de2 = fem::dsqrt(fem::pow2(ptwo(2, 5)) + fem::pow2(dpcm));
  //C
  de0 = fem::dsqrt(fem::pow2(px0) + fem::pow2(py0) + fem::pow2(pz0) +
    fem::pow2(xm0));
  dbex = px0 / de0;
  dbey = py0 / de0;
  dbez = pz0 / de0;
  //C     boost the reference frame up by beta (pznew=gam(pz+beta e)):
  beta2 = fem::pow2(dbex) + fem::pow2(dbey) + fem::pow2(dbez);
  gam = 1.e0 / fem::dsqrt(1.e0 - beta2);
  if (beta2 >= 0.9999999999999e0) {
    write(6, star), "1", dbex, dbey, dbez, beta2, gam;
  }
  //C
  lorenz(de1, dpx, dpy, dpz, -dbex, -dbey, -dbez);
  ptwo(1, 1) = pxnew;
  ptwo(1, 2) = pynew;
  ptwo(1, 3) = pznew;
  ptwo(1, 4) = enenew;
  lorenz(de2, -dpx, -dpy, -dpz, -dbex, -dbey, -dbez);
  ptwo(2, 1) = pxnew;
  ptwo(2, 2) = pynew;
  ptwo(2, 3) = pznew;
  ptwo(2, 4) = enenew;
  //C
}

struct htop_save
{
  float ftime;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int id;
  int idabs;
  int inozpc;
  int ipamax;
  int ipar;
  arr<int> it;
  int npar;
  double ptwox;
  double ptwoy;
  double ptwoz;
  float rnum;
  double xmdq;

  htop_save() :
    ftime(fem::float0),
    i(fem::int0),
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    i4(fem::int0),
    id(fem::int0),
    idabs(fem::int0),
    inozpc(fem::int0),
    ipamax(fem::int0),
    ipar(fem::int0),
    it(dimension(4), fem::fill0),
    npar(fem::int0),
    ptwox(fem::double0),
    ptwoy(fem::double0),
    ptwoz(fem::double0),
    rnum(fem::float0),
    xmdq(fem::double0)
  {}
};

//C
//C=======================================================================
void
htop(
  common& cmn)
{
  FEM_CMN_SVE(htop);
  common_write write(cmn);
  // COMMON hmain2
  const int maxstr = 150001;
  arr_cref<float, 2> patt(cmn.patt, dimension(maxstr, 4));
  // COMMON hmain1
  int& natt = cmn.natt;
  // COMMON prec1
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
  // COMMON ilist7
  arr_ref<int> lstrg0(cmn.lstrg0, dimension(maxptn));
  arr_ref<int> lpart0(cmn.lpart0, dimension(maxptn));
  // COMMON arprc
  arr_cref<int> itypar(cmn.itypar, dimension(maxstr));
  arr_cref<float> gxar(cmn.gxar, dimension(maxstr));
  arr_cref<float> gyar(cmn.gyar, dimension(maxstr));
  arr_cref<float> gzar(cmn.gzar, dimension(maxstr));
  arr_cref<float> ftar(cmn.ftar, dimension(maxstr));
  arr_cref<float> pxar(cmn.pxar, dimension(maxstr));
  arr_cref<float> pyar(cmn.pyar, dimension(maxstr));
  arr_cref<float> pzar(cmn.pzar, dimension(maxstr));
  arr_cref<float> pear(cmn.pear, dimension(maxstr));
  arr_cref<float> xmar(cmn.xmar, dimension(maxstr));
  // COMMON decom
  arr_ref<double, 2> ptwo(cmn.ptwo, dimension(2, 5));
  // COMMON noprec
  int& nnozpc = cmn.nnozpc;
  const int maxidl = 4001;
  arr_ref<int> itypn(cmn.itypn, dimension(maxidl));
  arr_ref<float> gxn(cmn.gxn, dimension(maxidl));
  arr_ref<float> gyn(cmn.gyn, dimension(maxidl));
  arr_ref<float> gzn(cmn.gzn, dimension(maxidl));
  arr_ref<float> ftn(cmn.ftn, dimension(maxidl));
  arr_ref<float> pxn(cmn.pxn, dimension(maxidl));
  arr_ref<float> pyn(cmn.pyn, dimension(maxidl));
  arr_ref<float> pzn(cmn.pzn, dimension(maxidl));
  arr_ref<float> een(cmn.een, dimension(maxidl));
  arr_ref<float> xmn(cmn.xmn, dimension(maxidl));
  // COMMON soft
  arr_ref<int> njsgs(cmn.njsgs, dimension(maxstr));
  // COMMON anim
  int& isoft = cmn.isoft;
  // COMMON precpa
  arr_ref<double> vxp0(cmn.vxp0, dimension(maxptn));
  arr_ref<double> vyp0(cmn.vyp0, dimension(maxptn));
  arr_ref<double> vzp0(cmn.vzp0, dimension(maxptn));
  arr_cref<double> xstrg0(cmn.xstrg0, dimension(maxptn));
  arr_cref<double> ystrg0(cmn.ystrg0, dimension(maxptn));
  arr_ref<double> xstrg(cmn.xstrg, dimension(maxptn));
  arr_ref<double> ystrg(cmn.ystrg, dimension(maxptn));
  arr_cref<int> istrg0(cmn.istrg0, dimension(maxptn));
  arr_ref<int> istrg(cmn.istrg, dimension(maxptn));
  // COMMON para7
  int& ioscar = cmn.ioscar;
  int& nsmbbbar = cmn.nsmbbbar;
  int& nsmmeson = cmn.nsmmeson;
  // COMMON snn
  float& epsipz = cmn.epsipz;
  float& epsipt = cmn.epsipt;
  float& pzproj = cmn.pzproj;
  float& pztarg = cmn.pztarg;
  //
  // SAVE
  float& ftime = sve.ftime;
  int& i = sve.i;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& i4 = sve.i4;
  int& id = sve.id;
  int& idabs = sve.idabs;
  int& inozpc = sve.inozpc;
  int& ipamax = sve.ipamax;
  int& ipar = sve.ipar;
  arr_ref<int> it(sve.it, dimension(4));
  int& npar = sve.npar;
  double& ptwox = sve.ptwox;
  double& ptwoy = sve.ptwoy;
  double& ptwoz = sve.ptwoz;
  float& rnum = sve.rnum;
  double& xmdq = sve.xmdq;
  //
  //C
  //Cc      SAVE /HMAIN2/
  //Cc      SAVE /HMAIN1/
  //Cc      SAVE /PARA1/
  //Cc      SAVE /prec1/
  //Cc      SAVE /ilist7/
  //Cc      SAVE /ARPRC/
  //Cc      SAVE /decom/
  //Cc      SAVE /RNDF77/
  //Cc      SAVE /NOPREC/
  //Cc      SAVE /HPARNT/
  //C     7/20/01: use double precision
  //C     otherwise sometimes beta>1 and gamma diverge in lorenz():
  //Cc      SAVE /SOFT/
  //Cc      SAVE /anim/
  //Clin-8/2015:
  //C      DOUBLE PRECISION  vxp0,vyp0,vzp0
  //C      common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
  //Cc      SAVE /precpa/
  //C
  npar = 0;
  nnozpc = 0;
  //Clin-5b/2008 calculate the number of hadrons to be converted to q/qbar:
  if ((isoft == 4 || isoft == 5) && (ioscar == 2 || ioscar == 3)) {
    nsmbbbar = 0;
    nsmmeson = 0;
    FEM_DO_SAFE(i, 1, natt) {
      id = itypar(i);
      idabs = fem::iabs(id);
      i2 = fem::mod(idabs / 10, 10);
      //Clin-9/2011 determine spectator nucleons consistently
      //C              if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
      //C     1             .ge.(HINT1(1)/2*0.99).and.
      //C     2             .and.(id.eq.2112.or.id.eq.2212)) then
      if (fem::abs(pxar(i)) <= epsipt && fem::abs(pyar(
          i)) <= epsipt && (pzar(i) > fem::amax1(0.f, pzproj -
          epsipz) || pzar(i) < (-pztarg + epsipz)) && (id == 2112 ||
          id == 2212)) {
        //C     spectator proj or targ nucleons without interactions, do not enter ZPC:
      }
      else if (idabs > 1000 && i2 != 0) {
        //C     baryons to be converted to q/qbar:
        nsmbbbar++;
      }
      else if ((idabs > 100 && idabs < 1000) || idabs > 10000) {
        //C     mesons to be converted to q/qbar:
        nsmmeson++;
      }
    }
    //C
    //Clin-6/2009:
    if (ioscar == 2 || ioscar == 3) {
      write(92, star), cmn.iaevt, cmn.miss, 3 * nsmbbbar + 2 * nsmmeson,
        nsmbbbar, nsmmeson, natt, natt - nsmbbbar - nsmmeson;
    }
    //C           write(92,*) iaevt, 3*nsmbbbar+2*nsmmeson
    //C           write(92,*) ' event#, total # of initial partons after string
    //C     1 melting'
    //C           write(92,*) 'String melting converts ',nsmbbbar, ' baryons &'
    //C     1, nsmmeson, ' mesons'
    //C           write(92,*) 'Total # of initial particles= ',natt
    //C           write(92,*) 'Total # of initial particles (gamma,e,muon,...)
    //C     1 not entering ZPC= ',natt-nsmbbbar-nsmmeson
  }
  //Clin-5b/2008-over
  FEM_DO_SAFE(i, 1, natt) {
    id = itypar(i);
    idabs = fem::iabs(id);
    i4 = fem::mod(idabs / 1000, 10);
    i3 = fem::mod(idabs / 100, 10);
    i2 = fem::mod(idabs / 10, 10);
    i1 = fem::mod(idabs, 10);
    rnum = ranart(cmn.nseed);
    ftime = 0.197f * pear(i) / (fem::pow2(pxar(i)) + fem::pow2(pyar(
      i)) + fem::pow2(xmar(i)));
    inozpc = 0;
    it(1) = 0;
    it(2) = 0;
    it(3) = 0;
    it(4) = 0;
    //C
    //Clin-9/2011 determine spectator nucleons consistently
    //C           if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
    //C     1 .ge.(HINT1(1)/2*0.99).and.((id.eq.2112).or.(id.eq.2212))) then
    if (fem::abs(pxar(i)) <= epsipt && fem::abs(pyar(i)) <= epsipt &&
        (pzar(i) > fem::amax1(0.f, pzproj - epsipz) || pzar(i) < (
        -pztarg + epsipz)) && (id == 2112 || id == 2212)) {
      //C     spectator proj or targ nucleons without interactions, do not enter ZPC:
      inozpc = 1;
    }
    else if (idabs > 1000 && i2 != 0) {
      //C     baryons:
      if (((i4 == 1 || i4 == 2) && i4 == i3) || (i4 == 3 && i3 == 3)) {
        if (i1 == 2) {
          if (rnum <= (1.f / 2.f)) {
            it(1) = i4;
            it(2) = i3 * 1000 + i2 * 100 + 1;
          }
          else if (rnum <= (2.f / 3.f)) {
            it(1) = i4;
            it(2) = i3 * 1000 + i2 * 100 + 3;
          }
          else {
            it(1) = i2;
            it(2) = i4 * 1000 + i3 * 100 + 3;
          }
        }
        else if (i1 == 4) {
          if (rnum <= (2.f / 3.f)) {
            it(1) = i4;
            it(2) = i3 * 1000 + i2 * 100 + 3;
          }
          else {
            it(1) = i2;
            it(2) = i4 * 1000 + i3 * 100 + 3;
          }
        }
      }
      else if (i4 == 1 || i4 == 2) {
        if (i1 == 2) {
          if (rnum <= (1.f / 2.f)) {
            it(1) = i2;
            it(2) = i4 * 1000 + i3 * 100 + 1;
          }
          else if (rnum <= (2.f / 3.f)) {
            it(1) = i2;
            it(2) = i4 * 1000 + i3 * 100 + 3;
          }
          else {
            it(1) = i4;
            it(2) = i3 * 1000 + i2 * 100 + 3;
          }
        }
        else if (i1 == 4) {
          if (rnum <= (2.f / 3.f)) {
            it(1) = i2;
            it(2) = i4 * 1000 + i3 * 100 + 3;
          }
          else {
            it(1) = i4;
            it(2) = i3 * 1000 + i2 * 100 + 3;
          }
        }
      }
      else if (i4 >= 3) {
        it(1) = i4;
        if (i3 < i2) {
          it(2) = i2 * 1000 + i3 * 100 + 1;
        }
        else {
          it(2) = i3 * 1000 + i2 * 100 + 3;
        }
      }
      //C       antibaryons:
      if (id < 0) {
        it(1) = -it(1);
        it(2) = -it(2);
      }
      //C     isoft=4or5 decompose diquark flavor it(2) to two quarks it(3)&(4):
      if (isoft == 4 || isoft == 5) {
        it(3) = fem::mod(it(2) / 1000, 10);
        it(4) = fem::mod(it(2) / 100, 10);
      }
      //C
    }
    else if ((idabs > 100 && idabs < 1000) || idabs > 10000) {
      //C     mesons:
      if (i3 == i2) {
        if (i3 == 1 || i3 == 2) {
          if (rnum <= 0.5f) {
            it(1) = 1;
            it(2) = -1;
          }
          else {
            it(1) = 2;
            it(2) = -2;
          }
        }
        else {
          it(1) = i3;
          it(2) = -i3;
        }
      }
      else {
        if ((fem::isign(1, id) * fem::pow((-1), i3)) == 1) {
          it(1) = i3;
          it(2) = -i2;
        }
        else {
          it(1) = i2;
          it(2) = -i3;
        }
      }
    }
    else {
      //C     save other particles (leptons and photons) outside of ZPC:
      inozpc = 1;
    }
    //C
    if (inozpc == 1) {
      njsgs(i) = 0;
      nnozpc++;
      itypn(nnozpc) = itypar(i);
      pxn(nnozpc) = pxar(i);
      pyn(nnozpc) = pyar(i);
      pzn(nnozpc) = pzar(i);
      een(nnozpc) = pear(i);
      xmn(nnozpc) = xmar(i);
      gxn(nnozpc) = gxar(i);
      gyn(nnozpc) = gyar(i);
      gzn(nnozpc) = gzar(i);
      ftn(nnozpc) = ftar(i);
    }
    else {
      njsgs(i) = 2;
      ptwo(1, 5) = fem::dble(ulmass(it(1)));
      ptwo(2, 5) = fem::dble(ulmass(it(2)));
      decomp(cmn, fem::dble(patt(i, 1)), fem::dble(patt(i, 2)),
        fem::dble(patt(i, 3)), fem::dble(xmar(i)), i, it(1));
      ipamax = 2;
      if ((isoft == 4 || isoft == 5) && fem::iabs(it(2)) > 1000) {
        ipamax = 1;
      }
      FEM_DO_SAFE(ipar, 1, ipamax) {
        npar++;
        ityp0(npar) = it(ipar);
        px0(npar) = ptwo(ipar, 1);
        py0(npar) = ptwo(ipar, 2);
        pz0(npar) = ptwo(ipar, 3);
        e0(npar) = ptwo(ipar, 4);
        xmass0(npar) = ptwo(ipar, 5);
        gx0(npar) = fem::dble(gxar(i));
        gy0(npar) = fem::dble(gyar(i));
        gz0(npar) = fem::dble(gzar(i));
        ft0(npar) = fem::dble(ftime);
        lstrg0(npar) = i;
        lpart0(npar) = ipar;
        vxp0(npar) = fem::dble(patt(i, 1) / patt(i, 4));
        vyp0(npar) = fem::dble(patt(i, 2) / patt(i, 4));
        vzp0(npar) = fem::dble(patt(i, 3) / patt(i, 4));
        //Clin-8/2015: set parent string information for this parton:
        xstrg(npar) = xstrg0(i);
        ystrg(npar) = ystrg0(i);
        istrg(npar) = istrg0(i);
      }
      //C
      if ((isoft == 4 || isoft == 5) && fem::iabs(it(2)) > 1000) {
        njsgs(i) = 3;
        xmdq = ptwo(2, 5);
        ptwo(1, 5) = fem::dble(ulmass(it(3)));
        ptwo(2, 5) = fem::dble(ulmass(it(4)));
        //C     8/19/02 avoid actual argument in common blocks of DECOMP:
        //C                 call decomp(ptwo(2,1),ptwo(2,2),ptwo(2,3),xmdq)
        ptwox = ptwo(2, 1);
        ptwoy = ptwo(2, 2);
        ptwoz = ptwo(2, 3);
        decomp(cmn, ptwox, ptwoy, ptwoz, xmdq, i, it(1));
        //C
        FEM_DO_SAFE(ipar, 1, 2) {
          npar++;
          ityp0(npar) = it(ipar + 2);
          px0(npar) = ptwo(ipar, 1);
          py0(npar) = ptwo(ipar, 2);
          pz0(npar) = ptwo(ipar, 3);
          e0(npar) = ptwo(ipar, 4);
          xmass0(npar) = ptwo(ipar, 5);
          gx0(npar) = fem::dble(gxar(i));
          gy0(npar) = fem::dble(gyar(i));
          gz0(npar) = fem::dble(gzar(i));
          ft0(npar) = fem::dble(ftime);
          lstrg0(npar) = i;
          lpart0(npar) = ipar + 1;
          vxp0(npar) = fem::dble(patt(i, 1) / patt(i, 4));
          vyp0(npar) = fem::dble(patt(i, 2) / patt(i, 4));
          vzp0(npar) = fem::dble(patt(i, 3) / patt(i, 4));
          //Clin-8/2015: set parent string information for this parton:
          xstrg(npar) = xstrg0(i);
          ystrg(npar) = ystrg0(i);
          istrg(npar) = istrg0(i);
        }
      }
      //C
    }
  }
  cmn.mul = npar;
  //C
  //Clin-5b/2008:
  if ((isoft == 4 || isoft == 5) && (ioscar == 2 || ioscar == 3)) {
    if ((natt - nsmbbbar - nsmmeson) != nnozpc) {
      write(92, star),
        "Problem with the total # of initial particles (gamma,e,muon,...) not "
        "entering ZPC";
    }
    if ((3 * nsmbbbar + 2 * nsmmeson) != npar) {
      write(92, star),
        "Problem with the total # of initial partons after string melting";
    }
  }
  //C
}

struct resmass_save
{
  float amass;
  float dm;
  float dmax;
  float dmin;
  float fm;
  float fmass;
  int ntry1;
  float wid;

  resmass_save() :
    amass(fem::float0),
    dm(fem::float0),
    dmax(fem::float0),
    dmin(fem::float0),
    fm(fem::float0),
    fmass(fem::float0),
    ntry1(fem::int0),
    wid(fem::float0)
  {}
};

//C
//C=======================================================================
//Clin-5/2011-add finite width to resonances (rho,omega,eta,K*,phi,Delta) after formation:
float
resmass(
  common& cmn,
  int const& kf)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(resmass);
  int& nseed = cmn.nseed;
  //
  float& amass = sve.amass;
  float& dm = sve.dm;
  float& dmax = sve.dmax;
  float& dmin = sve.dmin;
  float& fm = sve.fm;
  float& fmass = sve.fmass;
  int& ntry1 = sve.ntry1;
  float& wid = sve.wid;
  const float arho = 0.775f;
  const float wrho = 0.149f;
  const float aeta = 0.548f;
  const float weta = 1.30e-6f;
  const float aomega = 0.783f;
  const float womega = 0.00849f;
  const float aks = 0.894f;
  const float wks = 0.0498f;
  const float aphi = 1.019f;
  const float wphi = 0.00426f;
  const float adelta = 1.232f;
  const float wdelta = 0.118f;
  //C
  if (kf == 113 || fem::abs(kf) == 213) {
    amass = arho;
    wid = wrho;
  }
  else if (kf == 221) {
    amass = aeta;
    wid = weta;
  }
  else if (kf == 223) {
    amass = aomega;
    wid = womega;
  }
  else if (fem::abs(kf) == 313 || fem::abs(kf) == 323) {
    amass = aks;
    wid = wks;
  }
  else if (kf == 333) {
    amass = aphi;
    wid = wphi;
  }
  else if (fem::abs(kf) == 1114 || fem::abs(kf) == 2114 || fem::abs(
    kf) == 2214 || fem::abs(kf) == 2224) {
    amass = adelta;
    wid = wdelta;
  }
  dmin = amass - 2 * wid;
  dmax = amass + 2 * wid;
  //C     Delta mass needs to be big enough to decay to N+pi:
  if (amass == adelta) {
    dmin = 1.078f;
  }
  //C
  fm = 1.f;
  ntry1 = 0;
  statement_10:
  dm = ranart(nseed) * (dmax - dmin) + dmin;
  ntry1++;
  fmass = fem::pow2((amass * wid)) / (fem::pow2((fem::pow2(dm) -
    fem::pow2(amass))) + fem::pow2((amass * wid)));
  //Check      write (99,*) ntry1,kf,amass,wid,fmass,DM
  if ((ranart(nseed) > fmass / fm) && (ntry1 <= 10)) {
    goto statement_10;
  }
  //C
  return_value = dm;
  //C
  return return_value;
}

struct exchge_save
{
  double ft;
  double gx;
  double gy;
  double gz;
  int k1;
  int k2;
  double pe;
  double pm;
  double px;
  double py;
  double pz;

  exchge_save() :
    ft(fem::double0),
    gx(fem::double0),
    gy(fem::double0),
    gz(fem::double0),
    k1(fem::int0),
    k2(fem::int0),
    pe(fem::double0),
    pm(fem::double0),
    px(fem::double0),
    py(fem::double0),
    pz(fem::double0)
  {}
};

//C
//C=======================================================================
void
exchge(
  common& cmn,
  int const& isg,
  int const& ipi,
  int const& jsg,
  int const& ipj)
{
  FEM_CMN_SVE(exchge);
  // COMMON soft
  const int maxstr = 150001;
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
  //
  // SAVE
  double& ft = sve.ft;
  double& gx = sve.gx;
  double& gy = sve.gy;
  double& gz = sve.gz;
  int& k1 = sve.k1;
  int& k2 = sve.k2;
  double& pe = sve.pe;
  double& pm = sve.pm;
  double& px = sve.px;
  double& py = sve.py;
  double& pz = sve.pz;
  //
  //C
  //Cc      SAVE /SOFT/
  //C
  k1 = k1sgs(isg, ipi);
  k2 = k2sgs(isg, ipi);
  px = pxsgs(isg, ipi);
  py = pysgs(isg, ipi);
  pz = pzsgs(isg, ipi);
  pe = pesgs(isg, ipi);
  pm = pmsgs(isg, ipi);
  gx = gxsgs(isg, ipi);
  gy = gysgs(isg, ipi);
  gz = gzsgs(isg, ipi);
  ft = ftsgs(isg, ipi);
  k1sgs(isg, ipi) = k1sgs(jsg, ipj);
  k2sgs(isg, ipi) = k2sgs(jsg, ipj);
  pxsgs(isg, ipi) = pxsgs(jsg, ipj);
  pysgs(isg, ipi) = pysgs(jsg, ipj);
  pzsgs(isg, ipi) = pzsgs(jsg, ipj);
  pesgs(isg, ipi) = pesgs(jsg, ipj);
  pmsgs(isg, ipi) = pmsgs(jsg, ipj);
  gxsgs(isg, ipi) = gxsgs(jsg, ipj);
  gysgs(isg, ipi) = gysgs(jsg, ipj);
  gzsgs(isg, ipi) = gzsgs(jsg, ipj);
  ftsgs(isg, ipi) = ftsgs(jsg, ipj);
  k1sgs(jsg, ipj) = k1;
  k2sgs(jsg, ipj) = k2;
  pxsgs(jsg, ipj) = px;
  pysgs(jsg, ipj) = py;
  pzsgs(jsg, ipj) = pz;
  pesgs(jsg, ipj) = pe;
  pmsgs(jsg, ipj) = pm;
  gxsgs(jsg, ipj) = gx;
  gysgs(jsg, ipj) = gy;
  gzsgs(jsg, ipj) = gz;
  ftsgs(jsg, ipj) = ft;
  //C
}

struct locldr_save
{
  double beta2;
  double bex;
  double bey;
  double bez;
  double dt0;
  double etot;
  arr<double> ftp0;
  double gam;
  int iearly;
  int ilate;
  int imax;
  int imin;
  int istep;
  int j;
  arr<double> pep0;
  arr<double> pxp0;
  arr<double> pyp0;
  arr<double> pzp0;

  locldr_save() :
    beta2(fem::double0),
    bex(fem::double0),
    bey(fem::double0),
    bez(fem::double0),
    dt0(fem::double0),
    etot(fem::double0),
    ftp0(dimension(3), fem::fill0),
    gam(fem::double0),
    iearly(fem::int0),
    ilate(fem::int0),
    imax(fem::int0),
    imin(fem::int0),
    istep(fem::int0),
    j(fem::int0),
    pep0(dimension(3), fem::fill0),
    pxp0(dimension(3), fem::fill0),
    pyp0(dimension(3), fem::fill0),
    pzp0(dimension(3), fem::fill0)
  {}
};

//C
//C=======================================================================
void
locldr(
  common& cmn,
  int const& icall,
  double& drlocl)
{
  FEM_CMN_SVE(locldr);
  common_write write(cmn);
  // COMMON loclco
  arr_cref<double> gxp(cmn.gxp, dimension(3));
  arr_cref<double> gyp(cmn.gyp, dimension(3));
  arr_cref<double> gzp(cmn.gzp, dimension(3));
  arr_cref<double> ftp(cmn.ftp, dimension(3));
  arr_cref<double> pxp(cmn.pxp, dimension(3));
  arr_cref<double> pyp(cmn.pyp, dimension(3));
  arr_cref<double> pzp(cmn.pzp, dimension(3));
  arr_cref<double> pep(cmn.pep, dimension(3));
  arr_cref<double> pmp(cmn.pmp, dimension(3));
  // COMMON prtn23
  arr_ref<double> gxp0(cmn.gxp0, dimension(3));
  arr_ref<double> gyp0(cmn.gyp0, dimension(3));
  arr_ref<double> gzp0(cmn.gzp0, dimension(3));
  double& ft0fom = cmn.ft0fom;
  // COMMON lor
  double& enenew = cmn.enenew;
  double& pxnew = cmn.pxnew;
  double& pynew = cmn.pynew;
  double& pznew = cmn.pznew;
  //
  // SAVE
  double& beta2 = sve.beta2;
  double& bex = sve.bex;
  double& bey = sve.bey;
  double& bez = sve.bez;
  double& dt0 = sve.dt0;
  double& etot = sve.etot;
  arr_ref<double> ftp0(sve.ftp0, dimension(3));
  double& gam = sve.gam;
  int& iearly = sve.iearly;
  int& ilate = sve.ilate;
  int& imax = sve.imax;
  int& imin = sve.imin;
  int& istep = sve.istep;
  int& j = sve.j;
  arr_ref<double> pep0(sve.pep0, dimension(3));
  arr_ref<double> pxp0(sve.pxp0, dimension(3));
  arr_ref<double> pyp0(sve.pyp0, dimension(3));
  arr_ref<double> pzp0(sve.pzp0, dimension(3));
  //
  //C
  //Cc      SAVE /loclco/
  //Cc      SAVE /prtn23/
  //Cc      SAVE /lor/
  //C     for 2-body kinematics:
  if (icall == 2) {
    etot = pep(1) + pep(2);
    bex = (pxp(1) + pxp(2)) / etot;
    bey = (pyp(1) + pyp(2)) / etot;
    bez = (pzp(1) + pzp(2)) / etot;
    //C     boost the reference frame down by beta to get to the pair rest frame:
    FEM_DO_SAFE(j, 1, 2) {
      beta2 = fem::pow2(bex) + fem::pow2(bey) + fem::pow2(bez);
      gam = 1.e0 / fem::dsqrt(1.e0 - beta2);
      if (beta2 >= 0.9999999999999e0) {
        write(6, star), "4", pxp(1), pxp(2), pyp(1), pyp(2), pzp(1),
          pzp(2), pep(1), pep(2), pmp(1), pmp(2), fem::dsqrt(
          fem::pow2(pxp(1)) + fem::pow2(pyp(1)) + fem::pow2(pzp(1)) +
          fem::pow2(pmp(1))) / pep(1), fem::dsqrt(fem::pow2(pxp(1)) +
          fem::pow2(pyp(1)) + fem::pow2(pzp(1))) / pep(1);
        write(6, star), "4a", pxp(1) + pxp(2), pyp(1) + pyp(2), pzp(
          1) + pzp(2), etot;
        write(6, star), "4b", bex, bey, bez, beta2, gam;
      }
      //C
      lorenz(ftp(j), gxp(j), gyp(j), gzp(j), bex, bey, bez);
      gxp0(j) = pxnew;
      gyp0(j) = pynew;
      gzp0(j) = pznew;
      ftp0(j) = enenew;
      lorenz(pep(j), pxp(j), pyp(j), pzp(j), bex, bey, bez);
      pxp0(j) = pxnew;
      pyp0(j) = pynew;
      pzp0(j) = pznew;
      pep0(j) = enenew;
    }
    //C
    if (ftp0(1) >= ftp0(2)) {
      ilate = 1;
      iearly = 2;
    }
    else {
      ilate = 2;
      iearly = 1;
    }
    ft0fom = ftp0(ilate);
    //C
    dt0 = ftp0(ilate) - ftp0(iearly);
    gxp0(iearly) += pxp0(iearly) / pep0(iearly) * dt0;
    gyp0(iearly) += pyp0(iearly) / pep0(iearly) * dt0;
    gzp0(iearly) += pzp0(iearly) / pep0(iearly) * dt0;
    drlocl = fem::dsqrt(fem::pow2((gxp0(ilate) - gxp0(iearly))) +
      fem::pow2((gyp0(ilate) - gyp0(iearly))) + fem::pow2((gzp0(ilate) -
      gzp0(iearly))));
    //C     for 3-body kinematics, used for baryons formation:
  }
  else if (icall == 3) {
    etot = pep(1) + pep(2) + pep(3);
    bex = (pxp(1) + pxp(2) + pxp(3)) / etot;
    bey = (pyp(1) + pyp(2) + pyp(3)) / etot;
    bez = (pzp(1) + pzp(2) + pzp(3)) / etot;
    beta2 = fem::pow2(bex) + fem::pow2(bey) + fem::pow2(bez);
    gam = 1.e0 / fem::dsqrt(1.e0 - beta2);
    if (beta2 >= 0.9999999999999e0) {
      write(6, star), "5", bex, bey, bez, beta2, gam;
    }
    //C     boost the reference frame down by beta to get to the 3-parton rest frame:
    FEM_DO_SAFE(j, 1, 3) {
      lorenz(ftp(j), gxp(j), gyp(j), gzp(j), bex, bey, bez);
      gxp0(j) = pxnew;
      gyp0(j) = pynew;
      gzp0(j) = pznew;
      ftp0(j) = enenew;
      lorenz(pep(j), pxp(j), pyp(j), pzp(j), bex, bey, bez);
      pxp0(j) = pxnew;
      pyp0(j) = pynew;
      pzp0(j) = pznew;
      pep0(j) = enenew;
    }
    //C
    if (ftp0(1) > ftp0(2)) {
      ilate = 1;
      if (ftp0(3) > ftp0(1)) {
        ilate = 3;
      }
    }
    else {
      ilate = 2;
      if (ftp0(3) >= ftp0(2)) {
        ilate = 3;
      }
    }
    ft0fom = ftp0(ilate);
    //C
    if (ilate == 1) {
      imin = 2;
      imax = 3;
      istep = 1;
    }
    else if (ilate == 2) {
      imin = 1;
      imax = 3;
      istep = 2;
    }
    else if (ilate == 3) {
      imin = 1;
      imax = 2;
      istep = 1;
    }
    //C
    FEM_DOSTEP(iearly, imin, imax, istep) {
      dt0 = ftp0(ilate) - ftp0(iearly);
      gxp0(iearly) += pxp0(iearly) / pep0(iearly) * dt0;
      gyp0(iearly) += pyp0(iearly) / pep0(iearly) * dt0;
      gzp0(iearly) += pzp0(iearly) / pep0(iearly) * dt0;
    }
  }
  //C
}

struct coales_save
{
  static const int maxstr = 150001;

  double dp0;
  arr<double> dp1;
  double dplocl;
  double dr0;
  arr<double> dr1;
  double drlocl;
  int ibaryn;
  arr<int> iover;
  int ip;
  int ipi;
  int ipmax;
  int ipmin;
  int isg;
  int j;
  int jsg;

  coales_save() :
    dp0(fem::double0),
    dp1(dim1(2, 3), fem::fill0),
    dplocl(fem::double0),
    dr0(fem::double0),
    dr1(dim1(2, 3), fem::fill0),
    drlocl(fem::double0),
    ibaryn(fem::int0),
    iover(dimension(maxstr), fem::fill0),
    ip(fem::int0),
    ipi(fem::int0),
    ipmax(fem::int0),
    ipmin(fem::int0),
    isg(fem::int0),
    j(fem::int0),
    jsg(fem::int0)
  {}
};

const int coales_save::maxstr;

//C
//C=======================================================================
void
coales(
  common& cmn)
{
  FEM_CMN_SVE(coales);
  common_write write(cmn);
  const int maxstr = 150001;
  arr_cref<double, 2> pxsgs(cmn.pxsgs, dimension(maxstr, 3));
  arr_cref<double, 2> pysgs(cmn.pysgs, dimension(maxstr, 3));
  arr_cref<double, 2> pzsgs(cmn.pzsgs, dimension(maxstr, 3));
  arr_cref<double, 2> pesgs(cmn.pesgs, dimension(maxstr, 3));
  arr_cref<double, 2> pmsgs(cmn.pmsgs, dimension(maxstr, 3));
  arr_cref<double, 2> gxsgs(cmn.gxsgs, dimension(maxstr, 3));
  arr_cref<double, 2> gysgs(cmn.gysgs, dimension(maxstr, 3));
  arr_cref<double, 2> gzsgs(cmn.gzsgs, dimension(maxstr, 3));
  arr_cref<double, 2> ftsgs(cmn.ftsgs, dimension(maxstr, 3));
  arr_cref<int, 2> k2sgs(cmn.k2sgs, dimension(maxstr, 3));
  arr_cref<int> njsgs(cmn.njsgs, dimension(maxstr));
  double& dpcoal = cmn.dpcoal;
  double& drcoal = cmn.drcoal;
  arr_ref<double> gxp(cmn.gxp, dimension(3));
  arr_ref<double> gyp(cmn.gyp, dimension(3));
  arr_ref<double> gzp(cmn.gzp, dimension(3));
  arr_ref<double> ftp(cmn.ftp, dimension(3));
  arr_ref<double> pxp(cmn.pxp, dimension(3));
  arr_ref<double> pyp(cmn.pyp, dimension(3));
  arr_ref<double> pzp(cmn.pzp, dimension(3));
  arr_ref<double> pep(cmn.pep, dimension(3));
  arr_ref<double> pmp(cmn.pmp, dimension(3));
  int& nsg = cmn.nsg;
  //
  double& dp0 = sve.dp0;
  arr_ref<double> dp1(sve.dp1, dim1(2, 3));
  double& dplocl = sve.dplocl;
  double& dr0 = sve.dr0;
  arr_ref<double> dr1(sve.dr1, dim1(2, 3));
  double& drlocl = sve.drlocl;
  int& ibaryn = sve.ibaryn;
  arr_ref<int> iover(sve.iover, dimension(maxstr));
  int& ip = sve.ip;
  int& ipi = sve.ipi;
  int& ipmax = sve.ipmax;
  int& ipmin = sve.ipmin;
  int& isg = sve.isg;
  int& j = sve.j;
  int& jsg = sve.jsg;
  //C
  //Cc      SAVE /SOFT/
  //Cc      SAVE /coal/
  //Cc      SAVE /loclco/
  //Cc      SAVE /HJJET2/
  //C
  FEM_DO_SAFE(isg, 1, nsg) {
    iover(isg) = 0;
  }
  //C1     meson q coalesce with all available qbar:
  FEM_DO_SAFE(isg, 1, nsg) {
    if (njsgs(isg) != 2 || iover(isg) == 1) {
      goto statement_150;
    }
    //C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
    if (k2sgs(isg, 1) < 0) {
      write(6, star), "Antiquark appears in quark loop; stop";
      FEM_STOP(0);
    }
    //C
    FEM_DO_SAFE(j, 1, 2) {
      ftp(j) = ftsgs(isg, j);
      gxp(j) = gxsgs(isg, j);
      gyp(j) = gysgs(isg, j);
      gzp(j) = gzsgs(isg, j);
      pxp(j) = pxsgs(isg, j);
      pyp(j) = pysgs(isg, j);
      pzp(j) = pzsgs(isg, j);
      pmp(j) = pmsgs(isg, j);
      pep(j) = pesgs(isg, j);
    }
    locldr(cmn, 2, drlocl);
    dr0 = drlocl;
    //C     dp0^2 defined as (p1+p2)^2-(m1+m2)^2:
    dp0 = fem::dsqrt(2 * (pep(1) * pep(2) - pxp(1) * pxp(2) - pyp(
      1) * pyp(2) - pzp(1) * pzp(2) - pmp(1) * pmp(2)));
    //C
    FEM_DO_SAFE(jsg, 1, nsg) {
      //C     skip default or unavailable antiquarks:
      if (jsg == isg || iover(jsg) == 1) {
        goto statement_120;
      }
      if (njsgs(jsg) == 2) {
        ipmin = 2;
        ipmax = 2;
      }
      else if (njsgs(jsg) == 3 && k2sgs(jsg, 1) < 0) {
        ipmin = 1;
        ipmax = 3;
      }
      else {
        goto statement_120;
      }
      FEM_DO_SAFE(ip, ipmin, ipmax) {
        dplocl = fem::dsqrt(2 * (pep(1) * pesgs(jsg, ip) - pxp(1) * pxsgs(jsg,
          ip) - pyp(1) * pysgs(jsg, ip) - pzp(1) * pzsgs(jsg, ip) -
          pmp(1) * pmsgs(jsg, ip)));
        //C     skip if outside of momentum radius:
        if (dplocl > dpcoal) {
          goto statement_120;
        }
        ftp(2) = ftsgs(jsg, ip);
        gxp(2) = gxsgs(jsg, ip);
        gyp(2) = gysgs(jsg, ip);
        gzp(2) = gzsgs(jsg, ip);
        pxp(2) = pxsgs(jsg, ip);
        pyp(2) = pysgs(jsg, ip);
        pzp(2) = pzsgs(jsg, ip);
        pmp(2) = pmsgs(jsg, ip);
        pep(2) = pesgs(jsg, ip);
        locldr(cmn, 2, drlocl);
        //C     skip if outside of spatial radius:
        if (drlocl > drcoal) {
          goto statement_120;
        }
        //C     q_isg coalesces with qbar_jsg:
        if ((dp0 > dpcoal || dr0 > drcoal) || (drlocl < dr0)) {
          dp0 = dplocl;
          dr0 = drlocl;
          exchge(cmn, isg, 2, jsg, ip);
        }
      }
      statement_120:;
    }
    if (dp0 <= dpcoal && dr0 <= drcoal) {
      iover(isg) = 1;
    }
    statement_150:;
  }
  //C
  //C2     meson qbar coalesce with all available q:
  FEM_DO_SAFE(isg, 1, nsg) {
    if (njsgs(isg) != 2 || iover(isg) == 1) {
      goto statement_250;
    }
    //C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
    FEM_DO_SAFE(j, 1, 2) {
      ftp(j) = ftsgs(isg, j);
      gxp(j) = gxsgs(isg, j);
      gyp(j) = gysgs(isg, j);
      gzp(j) = gzsgs(isg, j);
      pxp(j) = pxsgs(isg, j);
      pyp(j) = pysgs(isg, j);
      pzp(j) = pzsgs(isg, j);
      pmp(j) = pmsgs(isg, j);
      pep(j) = pesgs(isg, j);
    }
    locldr(cmn, 2, drlocl);
    dr0 = drlocl;
    dp0 = fem::dsqrt(2 * (pep(1) * pep(2) - pxp(1) * pxp(2) - pyp(
      1) * pyp(2) - pzp(1) * pzp(2) - pmp(1) * pmp(2)));
    //C
    FEM_DO_SAFE(jsg, 1, nsg) {
      if (jsg == isg || iover(jsg) == 1) {
        goto statement_220;
      }
      if (njsgs(jsg) == 2) {
        ipmin = 1;
        ipmax = 1;
      }
      else if (njsgs(jsg) == 3 && k2sgs(jsg, 1) > 0) {
        ipmin = 1;
        ipmax = 3;
      }
      else {
        goto statement_220;
      }
      FEM_DO_SAFE(ip, ipmin, ipmax) {
        dplocl = fem::dsqrt(2 * (pep(2) * pesgs(jsg, ip) - pxp(2) * pxsgs(jsg,
          ip) - pyp(2) * pysgs(jsg, ip) - pzp(2) * pzsgs(jsg, ip) -
          pmp(2) * pmsgs(jsg, ip)));
        //C     skip if outside of momentum radius:
        if (dplocl > dpcoal) {
          goto statement_220;
        }
        ftp(1) = ftsgs(jsg, ip);
        gxp(1) = gxsgs(jsg, ip);
        gyp(1) = gysgs(jsg, ip);
        gzp(1) = gzsgs(jsg, ip);
        pxp(1) = pxsgs(jsg, ip);
        pyp(1) = pysgs(jsg, ip);
        pzp(1) = pzsgs(jsg, ip);
        pmp(1) = pmsgs(jsg, ip);
        pep(1) = pesgs(jsg, ip);
        locldr(cmn, 2, drlocl);
        //C     skip if outside of spatial radius:
        if (drlocl > drcoal) {
          goto statement_220;
        }
        //C     qbar_isg coalesces with q_jsg:
        if ((dp0 > dpcoal || dr0 > drcoal) || (drlocl < dr0)) {
          dp0 = dplocl;
          dr0 = drlocl;
          exchge(cmn, isg, 1, jsg, ip);
        }
      }
      statement_220:;
    }
    if (dp0 <= dpcoal && dr0 <= drcoal) {
      iover(isg) = 1;
    }
    statement_250:;
  }
  //C
  //C3     baryon q (antibaryon qbar) coalesce with all available q (qbar):
  FEM_DO_SAFE(isg, 1, nsg) {
    if (njsgs(isg) != 3 || iover(isg) == 1) {
      goto statement_350;
    }
    ibaryn = k2sgs(isg, 1);
    //C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
    FEM_DO_SAFE(j, 1, 2) {
      ftp(j) = ftsgs(isg, j);
      gxp(j) = gxsgs(isg, j);
      gyp(j) = gysgs(isg, j);
      gzp(j) = gzsgs(isg, j);
      pxp(j) = pxsgs(isg, j);
      pyp(j) = pysgs(isg, j);
      pzp(j) = pzsgs(isg, j);
      pmp(j) = pmsgs(isg, j);
      pep(j) = pesgs(isg, j);
    }
    locldr(cmn, 2, drlocl);
    dr1(2) = drlocl;
    dp1(2) = fem::dsqrt(2 * (pep(1) * pep(2) - pxp(1) * pxp(2) - pyp(
      1) * pyp(2) - pzp(1) * pzp(2) - pmp(1) * pmp(2)));
    //C
    ftp(2) = ftsgs(isg, 3);
    gxp(2) = gxsgs(isg, 3);
    gyp(2) = gysgs(isg, 3);
    gzp(2) = gzsgs(isg, 3);
    pxp(2) = pxsgs(isg, 3);
    pyp(2) = pysgs(isg, 3);
    pzp(2) = pzsgs(isg, 3);
    pmp(2) = pmsgs(isg, 3);
    pep(2) = pesgs(isg, 3);
    locldr(cmn, 2, drlocl);
    dr1(3) = drlocl;
    dp1(3) = fem::dsqrt(2 * (pep(1) * pep(2) - pxp(1) * pxp(2) - pyp(
      1) * pyp(2) - pzp(1) * pzp(2) - pmp(1) * pmp(2)));
    //C
    FEM_DO_SAFE(jsg, 1, nsg) {
      if (jsg == isg || iover(jsg) == 1) {
        goto statement_320;
      }
      if (njsgs(jsg) == 2) {
        if (ibaryn > 0) {
          ipmin = 1;
        }
        else {
          ipmin = 2;
        }
        ipmax = ipmin;
      }
      else if (njsgs(jsg) == 3 && (ibaryn * k2sgs(jsg, 1)) > 0) {
        ipmin = 1;
        ipmax = 3;
      }
      else {
        goto statement_320;
      }
      FEM_DO_SAFE(ip, ipmin, ipmax) {
        dplocl = fem::dsqrt(2 * (pep(1) * pesgs(jsg, ip) - pxp(1) * pxsgs(jsg,
          ip) - pyp(1) * pysgs(jsg, ip) - pzp(1) * pzsgs(jsg, ip) -
          pmp(1) * pmsgs(jsg, ip)));
        //C     skip if outside of momentum radius:
        if (dplocl > dpcoal) {
          goto statement_320;
        }
        ftp(2) = ftsgs(jsg, ip);
        gxp(2) = gxsgs(jsg, ip);
        gyp(2) = gysgs(jsg, ip);
        gzp(2) = gzsgs(jsg, ip);
        pxp(2) = pxsgs(jsg, ip);
        pyp(2) = pysgs(jsg, ip);
        pzp(2) = pzsgs(jsg, ip);
        pmp(2) = pmsgs(jsg, ip);
        pep(2) = pesgs(jsg, ip);
        locldr(cmn, 2, drlocl);
        //C     skip if outside of spatial radius:
        if (drlocl > drcoal) {
          goto statement_320;
        }
        //C     q_isg may coalesce with q_jsg for a baryon:
        ipi = 0;
        if (dp1(2) > dpcoal || dr1(2) > drcoal) {
          ipi = 2;
          if ((dp1(3) > dpcoal || dr1(3) > drcoal) && dr1(3) > dr1(2)) {
            ipi = 3;
          }
        }
        else if (dp1(3) > dpcoal || dr1(3) > drcoal) {
          ipi = 3;
        }
        else if (dr1(2) < dr1(3)) {
          if (drlocl < dr1(3)) {
            ipi = 3;
          }
        }
        else if (dr1(3) <= dr1(2)) {
          if (drlocl < dr1(2)) {
            ipi = 2;
          }
        }
        if (ipi != 0) {
          dp1(ipi) = dplocl;
          dr1(ipi) = drlocl;
          exchge(cmn, isg, ipi, jsg, ip);
        }
      }
      statement_320:;
    }
    if (dp1(2) <= dpcoal && dr1(2) <= drcoal && dp1(3) <= dpcoal &&
        dr1(3) <= drcoal) {
      iover(isg) = 1;
    }
    statement_350:;
  }
  //C
}

struct ptoh_save
{
  static const int maxstr = 150001;

  double beta2;
  double bex;
  double bey;
  double bez;
  double drlocl;
  double e1;
  double e2;
  double e3;
  double etot;
  double ftavg0;
  double gam;
  double gxavg0;
  double gyavg0;
  double gzavg0;
  int i;
  int ibs;
  int idqspn;
  int imspin;
  int inatt;
  arr<int> indx;
  int ipartn;
  int isg;
  int iuudd;
  int ix;
  int k1;
  int k1abs;
  int k2;
  int k2abs;
  int k3;
  int k3abs;
  int kdq;
  int kf;
  int kf1;
  int kf2;
  int ki;
  int kj;
  int kk;
  int kmax;
  int kmin;
  int ktemp;
  int mstj24;
  arr<int> ndiag;
  int npi0;
  int npich;
  int nrhoch;
  int nuudd;
  double p1;
  double p2;
  double p3;
  float ppi0;
  float prho0;
  double px1;
  double px2;
  double px3;
  double py1;
  double py2;
  double py3;
  double pz1;
  double pz2;
  double pz3;
  float tau0;
  arr<double> xmdiag;
  double xmpair;

  ptoh_save() :
    beta2(fem::double0),
    bex(fem::double0),
    bey(fem::double0),
    bez(fem::double0),
    drlocl(fem::double0),
    e1(fem::double0),
    e2(fem::double0),
    e3(fem::double0),
    etot(fem::double0),
    ftavg0(fem::double0),
    gam(fem::double0),
    gxavg0(fem::double0),
    gyavg0(fem::double0),
    gzavg0(fem::double0),
    i(fem::int0),
    ibs(fem::int0),
    idqspn(fem::int0),
    imspin(fem::int0),
    inatt(fem::int0),
    indx(dimension(maxstr), fem::fill0),
    ipartn(fem::int0),
    isg(fem::int0),
    iuudd(fem::int0),
    ix(fem::int0),
    k1(fem::int0),
    k1abs(fem::int0),
    k2(fem::int0),
    k2abs(fem::int0),
    k3(fem::int0),
    k3abs(fem::int0),
    kdq(fem::int0),
    kf(fem::int0),
    kf1(fem::int0),
    kf2(fem::int0),
    ki(fem::int0),
    kj(fem::int0),
    kk(fem::int0),
    kmax(fem::int0),
    kmin(fem::int0),
    ktemp(fem::int0),
    mstj24(fem::int0),
    ndiag(dimension(maxstr), fem::fill0),
    npi0(fem::int0),
    npich(fem::int0),
    nrhoch(fem::int0),
    nuudd(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    ppi0(fem::float0),
    prho0(fem::float0),
    px1(fem::double0),
    px2(fem::double0),
    px3(fem::double0),
    py1(fem::double0),
    py2(fem::double0),
    py3(fem::double0),
    pz1(fem::double0),
    pz2(fem::double0),
    pz3(fem::double0),
    tau0(fem::float0),
    xmdiag(dimension(maxstr), fem::fill0),
    xmpair(fem::double0)
  {}
};

const int ptoh_save::maxstr;

//C
//C=======================================================================
void
ptoh(
  common& cmn)
{
  FEM_CMN_SVE(ptoh);
  common_write write(cmn);
  // COMMON loclco
  arr_ref<double> gxp(cmn.gxp, dimension(3));
  arr_ref<double> gyp(cmn.gyp, dimension(3));
  arr_ref<double> gzp(cmn.gzp, dimension(3));
  arr_ref<double> ftp(cmn.ftp, dimension(3));
  arr_ref<double> pxp(cmn.pxp, dimension(3));
  arr_ref<double> pyp(cmn.pyp, dimension(3));
  arr_ref<double> pzp(cmn.pzp, dimension(3));
  arr_ref<double> pep(cmn.pep, dimension(3));
  arr_ref<double> pmp(cmn.pmp, dimension(3));
  // COMMON hmain1
  float& eatt = cmn.eatt;
  int& natt = cmn.natt;
  // COMMON hmain2
  const int maxstr = 150001;
  arr_ref<int, 2> katt(cmn.katt, dimension(maxstr, 4));
  arr_ref<float, 2> patt(cmn.patt, dimension(maxstr, 4));
  // COMMON hjjet2
  int& nsg = cmn.nsg;
  // COMMON arprnt
  arr_cref<float> arpar1(cmn.arpar1, dimension(100));
  // COMMON arprc
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
  // COMMON soft
  arr_cref<double, 2> pxsgs(cmn.pxsgs, dimension(maxstr, 3));
  arr_cref<double, 2> pysgs(cmn.pysgs, dimension(maxstr, 3));
  arr_cref<double, 2> pzsgs(cmn.pzsgs, dimension(maxstr, 3));
  arr_cref<double, 2> pesgs(cmn.pesgs, dimension(maxstr, 3));
  arr_cref<double, 2> pmsgs(cmn.pmsgs, dimension(maxstr, 3));
  arr_cref<double, 2> gxsgs(cmn.gxsgs, dimension(maxstr, 3));
  arr_cref<double, 2> gysgs(cmn.gysgs, dimension(maxstr, 3));
  arr_cref<double, 2> gzsgs(cmn.gzsgs, dimension(maxstr, 3));
  arr_cref<double, 2> ftsgs(cmn.ftsgs, dimension(maxstr, 3));
  arr_cref<int, 2> k2sgs(cmn.k2sgs, dimension(maxstr, 3));
  arr_cref<int> njsgs(cmn.njsgs, dimension(maxstr));
  // COMMON rndf77
  int& nseed = cmn.nseed;
  // COMMON anim
  int& isoft = cmn.isoft;
  // COMMON prtn23
  arr_cref<double> gxp0(cmn.gxp0, dimension(3));
  arr_cref<double> gyp0(cmn.gyp0, dimension(3));
  arr_cref<double> gzp0(cmn.gzp0, dimension(3));
  // COMMON ludat1
  arr_ref<int> mstj(cmn.mstj, dimension(200));
  // COMMON para7
  int& ioscar = cmn.ioscar;
  int& nsmbbbar = cmn.nsmbbbar;
  int& nsmmeson = cmn.nsmmeson;
  //
  // SAVE
  double& beta2 = sve.beta2;
  double& bex = sve.bex;
  double& bey = sve.bey;
  double& bez = sve.bez;
  double& e1 = sve.e1;
  double& e2 = sve.e2;
  double& e3 = sve.e3;
  double& etot = sve.etot;
  double& ftavg0 = sve.ftavg0;
  double& gam = sve.gam;
  double& gxavg0 = sve.gxavg0;
  double& gyavg0 = sve.gyavg0;
  double& gzavg0 = sve.gzavg0;
  int& i = sve.i;
  int& ibs = sve.ibs;
  int& idqspn = sve.idqspn;
  int& imspin = sve.imspin;
  int& inatt = sve.inatt;
  arr_ref<int> indx(sve.indx, dimension(maxstr));
  int& ipartn = sve.ipartn;
  int& isg = sve.isg;
  int& iuudd = sve.iuudd;
  int& ix = sve.ix;
  int& k1 = sve.k1;
  int& k1abs = sve.k1abs;
  int& k2 = sve.k2;
  int& k2abs = sve.k2abs;
  int& k3 = sve.k3;
  int& k3abs = sve.k3abs;
  int& kdq = sve.kdq;
  int& kf = sve.kf;
  int& kf1 = sve.kf1;
  int& kf2 = sve.kf2;
  int& ki = sve.ki;
  int& kj = sve.kj;
  int& kk = sve.kk;
  int& kmax = sve.kmax;
  int& kmin = sve.kmin;
  int& ktemp = sve.ktemp;
  int& mstj24 = sve.mstj24;
  arr_ref<int> ndiag(sve.ndiag, dimension(maxstr));
  int& npi0 = sve.npi0;
  int& npich = sve.npich;
  int& nrhoch = sve.nrhoch;
  int& nuudd = sve.nuudd;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  float& ppi0 = sve.ppi0;
  float& prho0 = sve.prho0;
  double& px1 = sve.px1;
  double& px2 = sve.px2;
  double& px3 = sve.px3;
  double& py1 = sve.py1;
  double& py2 = sve.py2;
  double& py3 = sve.py3;
  double& pz1 = sve.pz1;
  double& pz2 = sve.pz2;
  double& pz3 = sve.pz3;
  float& tau0 = sve.tau0;
  arr_ref<double> xmdiag(sve.xmdiag, dimension(maxstr));
  double& xmpair = sve.xmpair;
  //
  static const char* format_312 = "(i6,4(1x,f10.3),1x,i6,1x,i6)";
  //C
  //Clin-9/2012: improve precision for argument in sqrt():
  //Cc      SAVE /loclco/
  //Cc      SAVE /HMAIN1/
  //Cc      SAVE /HMAIN2/
  //Cc      SAVE /HJJET2/
  //Cc      SAVE /ARPRNT/
  //Cc      SAVE /ARPRC/
  //Cc      SAVE /SOFT/
  //Cc      SAVE /RNDF77/
  //Cc      SAVE /anim/
  //Cc      SAVE /prtn23/
  //Cc      SAVE /nzpc/
  //Cc      SAVE /lor/
  //Cc      SAVE /LUDAT1/
  //Clin 4/19/2006
  //Clin-5/2011
  //C
  coales(cmn);
  //C     obtain particle mass here without broadening by Breit-Wigner width:
  mstj24 = mstj(24);
  mstj(24) = 0;
  nuudd = 0;
  npich = 0;
  nrhoch = 0;
  ppi0 = 1.f;
  prho0 = 0.f;
  //C     determine hadron flavor except (pi0,rho0,eta,omega):
  FEM_DO_SAFE(isg, 1, nsg) {
    if (njsgs(isg) != 0) {
      natt++;
      k1 = k2sgs(isg, 1);
      k1abs = fem::iabs(k1);
      px1 = pxsgs(isg, 1);
      py1 = pysgs(isg, 1);
      pz1 = pzsgs(isg, 1);
      k2 = k2sgs(isg, 2);
      k2abs = fem::iabs(k2);
      px2 = pxsgs(isg, 2);
      py2 = pysgs(isg, 2);
      pz2 = pzsgs(isg, 2);
      //C     5/02/01 try lowest spin states as first choices,
      //C     i.e. octet baryons and pseudoscalar mesons (ibs=2*baryonspin+1):
      e1 = pesgs(isg, 1);
      e2 = pesgs(isg, 2);
      xmpair = fem::dsqrt(fem::pow2((e1 + e2)) - fem::pow2((px1 +
        px2)) - fem::pow2((py1 + py2)) - fem::pow2((pz1 + pz2)));
      ibs = 2;
      imspin = 0;
      if (k1 ==  - k2 && fem::iabs(k1) <= 2 && njsgs(isg) == 2) {
        nuudd++;
        xmdiag(nuudd) = xmpair;
        ndiag(nuudd) = natt;
      }
      k3 = 0;
      if ((isoft == 4 || isoft == 5) && njsgs(isg) == 3) {
        k3 = k2sgs(isg, 3);
        k3abs = fem::iabs(k3);
        px3 = pxsgs(isg, 3);
        py3 = pysgs(isg, 3);
        pz3 = pzsgs(isg, 3);
        e3 = pesgs(isg, 3);
        xmpair = fem::dsqrt(fem::pow2((e1 + e2 + e3)) - fem::pow2((
          px1 + px2 + px3)) - fem::pow2((py1 + py2 + py3)) -
          fem::pow2((pz1 + pz2 + pz3)));
      }
      //C*****     isoft=3 baryon decomposition is different:
      if (isoft == 3 && (k1abs > 1000 || k2abs > 1000)) {
        if (k1abs > 1000) {
          kdq = k1abs;
          kk = k2abs;
        }
        else {
          kdq = k2abs;
          kk = k1abs;
        }
        ki = fem::mod(kdq / 1000, 10);
        kj = fem::mod(kdq / 100, 10);
        if (fem::mod(kdq, 10) == 1) {
          idqspn = 0;
        }
        else {
          idqspn = 1;
        }
        //C
        if (kk > ki) {
          ktemp = kk;
          kk = kj;
          kj = ki;
          ki = ktemp;
        }
        else if (kk > kj) {
          ktemp = kk;
          kk = kj;
          kj = ktemp;
        }
        //C
        if (ki != kj && ki != kk && kj != kk) {
          if (idqspn == 0) {
            kf = 1000 * ki + 100 * kk + 10 * kj + ibs;
          }
          else {
            kf = 1000 * ki + 100 * kj + 10 * kk + ibs;
          }
        }
        else if (ki == kj && ki == kk) {
          //C     can only be decuplet baryons:
          kf = 1000 * ki + 100 * kj + 10 * kk + 4;
        }
        else {
          kf = 1000 * ki + 100 * kj + 10 * kk + ibs;
        }
        //C     form a decuplet baryon if the q+diquark mass is closer to its mass
        //C     (and if the diquark has spin 1):
        //Cc     for now only include Delta, which is present in ART:
        //Cc                 if(idqspn.eq.1.and.MOD(kf,10).eq.2) then
        if (kf == 2112 || kf == 2212) {
          if (fem::abs(fem::sngl(xmpair) - ulmass(kf)) > fem::abs(
              fem::sngl(xmpair) - ulmass(kf + 2))) {
            kf += 2;
          }
        }
        if (k1 < 0) {
          kf = -kf;
        }
        //Clin-6/22/01 isoft=4or5 baryons:
      }
      else if ((isoft == 4 || isoft == 5) && njsgs(isg) == 3) {
        if (k1abs > k2abs) {
          ki = k1abs;
          kk = k2abs;
        }
        else {
          ki = k2abs;
          kk = k1abs;
        }
        if (k3abs > ki) {
          kj = ki;
          ki = k3abs;
        }
        else if (k3abs < kk) {
          kj = kk;
          kk = k3abs;
        }
        else {
          kj = k3abs;
        }
        //C
        if (ki == kj && ki == kk) {
          //C     can only be decuplet baryons (Delta-,++, Omega):
          ibs = 4;
          kf = 1000 * ki + 100 * kj + 10 * kk + ibs;
        }
        else if (ki != kj && ki != kk && kj != kk) {
          //C     form Lambda or Sigma according to 3-quark mass,
          //C     for now neglect decuplet (Sigma*0 etc) which is absent in ART:
          ibs = 2;
          kf1 = 1000 * ki + 100 * kj + 10 * kk + ibs;
          kf2 = 1000 * ki + 100 * kk + 10 * kj + ibs;
          kf = kf1;
          if (fem::abs(fem::sngl(xmpair) - ulmass(kf1)) > fem::abs(
              fem::sngl(xmpair) - ulmass(kf2))) {
            kf = kf2;
          }
        }
        else {
          ibs = 2;
          kf = 1000 * ki + 100 * kj + 10 * kk + ibs;
          //Cc     for now only include Delta0,+ as decuplets, which are present in ART:
          if (kf == 2112 || kf == 2212) {
            if (fem::abs(fem::sngl(xmpair) - ulmass(kf)) > fem::abs(
                fem::sngl(xmpair) - ulmass(kf + 2))) {
              kf += 2;
            }
          }
        }
        if (k1 < 0) {
          kf = -kf;
        }
        //C*****     mesons:
      }
      else {
        if (k1abs == k2abs) {
          if (k1abs <= 2) {
            //C     treat diagonal mesons later in the subroutine:
            kf = 0;
          }
          else if (k1abs <= 3) {
            //C     do not form eta', only form phi from s-sbar, since no eta' in ART:
            kf = 333;
          }
          else {
            kf = 100 * k1abs + 10 * k1abs + 2 * imspin + 1;
          }
        }
        else {
          if (k1abs > k2abs) {
            kmax = k1abs;
            kmin = k2abs;
          }
          else if (k1abs < k2abs) {
            kmax = k2abs;
            kmin = k1abs;
          }
          kf = (100 * kmax + 10 * kmin + 2 * imspin + 1) * fem::isign(1,
            k1 + k2) * fem::pow((-1), kmax);
          //C     form a vector meson if the q+qbar mass is closer to its mass:
          if (fem::mod(fem::iabs(kf), 10) == 1) {
            if (fem::abs(fem::sngl(xmpair) - ulmass(fem::iabs(
                kf))) > fem::abs(fem::sngl(xmpair) - ulmass(fem::iabs(
                kf) + 2))) {
              kf = (fem::iabs(kf) + 2) * fem::isign(1, kf);
            }
          }
        }
      }
      itypar(natt) = kf;
      katt(natt, 1) = kf;
      if (fem::iabs(kf) == 211) {
        npich++;
      }
      else if (fem::iabs(kf) == 213) {
        nrhoch++;
      }
    }
    //Clin-7/2011-check charm hadron flavors:
    //C           if(k1abs.eq.4.or.k2abs.eq.4) then
    //C              if(k3.eq.0) then
    //C                 write(99,*) iaevt,k1,k2,kf,xmpair,
    //C     1                ULMASS(iabs(kf)),ULMASS(iabs(kf)+2),isg
    //C              else
    //C                 write(99,*) iaevt,k1,k2,k3,kf,xmpair,
    //C     1                ULMASS(iabs(kf)),ULMASS(iabs(kf)+2),isg
    //C              endif
    //C           endif
    //Clin-7/2011-end
  }
  //C     assume Npi0=(Npi+ + Npi-)/2, Nrho0=(Nrho+ + Nrho-)/2 on the average:
  if (nuudd != 0) {
    ppi0 = fem::ffloat(npich / 2) / fem::ffloat(nuudd);
    prho0 = fem::ffloat(nrhoch / 2) / fem::ffloat(nuudd);
  }
  //C     determine diagonal mesons (pi0,rho0,eta and omega) from uubar/ddbar:
  npi0 = 0;
  FEM_DO_SAFE(isg, 1, nsg) {
    if (k2sgs(isg, 1) ==  - k2sgs(isg, 2) && fem::iabs(k2sgs(isg,
        1)) <= 2 && njsgs(isg) == 2) {
      if (ranart(nseed) <= ppi0) {
        npi0++;
      }
    }
  }
  //C
  if (nuudd > 1) {
    index1(maxstr, nuudd, xmdiag, indx);
  }
  else {
    indx(1) = 1;
  }
  //C
  FEM_DO_SAFE(ix, 1, nuudd) {
    iuudd = indx(ix);
    inatt = ndiag(iuudd);
    if (ix <= npi0) {
      kf = 111;
    }
    else if (ranart(nseed) <= (prho0 / (1 - ppi0 + 0.00001f))) {
      kf = 113;
    }
    else {
      //C     at T=150MeV, thermal weights for eta and omega(spin1) are about the same:
      if (ranart(nseed) <= 0.5f) {
        kf = 221;
      }
      else {
        kf = 223;
      }
    }
    itypar(inatt) = kf;
    katt(inatt, 1) = kf;
  }
  //C  determine hadron formation time, position and momentum:
  inatt = 0;
  //Clin-6/2009 write out parton info after coalescence:
  if (ioscar == 3) {
    write(85, "(4i8,f10.4,5i5)"), cmn.iaevt, 3 * nsmbbbar + 2 * nsmmeson,
      nsmbbbar, nsmmeson, cmn.bimp, cmn.nelp, cmn.ninp, cmn.nelt,
      cmn.ninthj, cmn.miss;
  }
  //C
  FEM_DO_SAFE(isg, 1, nsg) {
    if (njsgs(isg) != 0) {
      inatt++;
      k1 = k2sgs(isg, 1);
      k1abs = fem::iabs(k1);
      px1 = pxsgs(isg, 1);
      py1 = pysgs(isg, 1);
      pz1 = pzsgs(isg, 1);
      k2 = k2sgs(isg, 2);
      k2abs = fem::iabs(k2);
      px2 = pxsgs(isg, 2);
      py2 = pysgs(isg, 2);
      pz2 = pzsgs(isg, 2);
      e1 = pesgs(isg, 1);
      e2 = pesgs(isg, 2);
      //C
      if (njsgs(isg) == 2) {
        pxar(inatt) = fem::sngl(px1 + px2);
        pyar(inatt) = fem::sngl(py1 + py2);
        pzar(inatt) = fem::sngl(pz1 + pz2);
        patt(inatt, 1) = pxar(inatt);
        patt(inatt, 2) = pyar(inatt);
        patt(inatt, 3) = pzar(inatt);
        etot = e1 + e2;
        //Clin-9/2012: improve precision for argument in sqrt():
        p1 = px1 + px2;
        p2 = py1 + py2;
        p3 = pz1 + pz2;
        //C
      }
      else if ((isoft == 4 || isoft == 5) && njsgs(isg) == 3) {
        px3 = pxsgs(isg, 3);
        py3 = pysgs(isg, 3);
        pz3 = pzsgs(isg, 3);
        e3 = pesgs(isg, 3);
        pxar(inatt) = fem::sngl(px1 + px2 + px3);
        pyar(inatt) = fem::sngl(py1 + py2 + py3);
        pzar(inatt) = fem::sngl(pz1 + pz2 + pz3);
        patt(inatt, 1) = pxar(inatt);
        patt(inatt, 2) = pyar(inatt);
        patt(inatt, 3) = pzar(inatt);
        etot = e1 + e2 + e3;
        //Clin-9/2012: improve precision for argument in sqrt():
        p1 = px1 + px2 + px3;
        p2 = py1 + py2 + py3;
        p3 = pz1 + pz2 + pz3;
        //C
      }
      xmar(inatt) = ulmass(itypar(inatt));
      //Clin-5/2011-add finite width to resonances (rho,omega,eta,K*,phi,Delta) after formation:
      kf = katt(inatt, 1);
      if (kf == 113 || fem::abs(kf) == 213 || kf == 221 ||
          kf == 223 || fem::abs(kf) == 313 || fem::abs(kf) == 323 ||
          kf == 333 || fem::abs(kf) == 1114 || fem::abs(
          kf) == 2114 || fem::abs(kf) == 2214 || fem::abs(
          kf) == 2224) {
        xmar(inatt) = resmass(cmn, kf);
      }
      //C
      pear(inatt) = fem::sqrt(fem::pow2(pxar(inatt)) + fem::pow2(pyar(
        inatt)) + fem::pow2(pzar(inatt)) + fem::pow2(xmar(inatt)));
      patt(inatt, 4) = pear(inatt);
      eatt += pear(inatt);
      ipartn = njsgs(isg);
      FEM_DO_SAFE(i, 1, ipartn) {
        ftp(i) = ftsgs(isg, i);
        gxp(i) = gxsgs(isg, i);
        gyp(i) = gysgs(isg, i);
        gzp(i) = gzsgs(isg, i);
        pxp(i) = pxsgs(isg, i);
        pyp(i) = pysgs(isg, i);
        pzp(i) = pzsgs(isg, i);
        pmp(i) = pmsgs(isg, i);
        pep(i) = pesgs(isg, i);
      }
      locldr(cmn, ipartn, sve.drlocl);
      //C
      tau0 = arpar1(1);
      ftavg0 = cmn.ft0fom + fem::dble(tau0);
      gxavg0 = 0e0;
      gyavg0 = 0e0;
      gzavg0 = 0e0;
      FEM_DO_SAFE(i, 1, ipartn) {
        gxavg0 += gxp0(i) / ipartn;
        gyavg0 += gyp0(i) / ipartn;
        gzavg0 += gzp0(i) / ipartn;
      }
      //Clin-9/2012: improve precision for argument in sqrt():
      //C            bex=dble(PXAR(inatt))/etot
      //C            bey=dble(PYAR(inatt))/etot
      //C            bez=dble(PZAR(inatt))/etot
      bex = p1 / etot;
      bey = p2 / etot;
      bez = p3 / etot;
      //C
      beta2 = fem::pow2(bex) + fem::pow2(bey) + fem::pow2(bez);
      gam = 1.e0 / fem::dsqrt(1.e0 - beta2);
      if (beta2 >= 0.9999999999999e0) {
        write(6, star), "2", bex, bey, bez, beta2, gam;
      }
      //C
      lorenz(ftavg0, gxavg0, gyavg0, gzavg0, -bex, -bey, -bez);
      gxar(inatt) = fem::sngl(cmn.pxnew);
      gyar(inatt) = fem::sngl(cmn.pynew);
      gzar(inatt) = fem::sngl(cmn.pznew);
      ftar(inatt) = fem::sngl(cmn.enenew);
      //Clin 4/19/2006 write out parton info after coalescence:
      if (ioscar == 3) {
        write(85, "(i6,4(1x,f10.3),1x,i6,1x,i6,1x,f10.3)"), k2sgs(isg,
          1), px1, py1, pz1, pmsgs(isg, 1), inatt, katt(inatt, 1),
          xmar(inatt);
        write(85, format_312), k2sgs(isg, 2), px2, py2, pz2, pmsgs(isg,
          2), inatt, katt(inatt, 1);
        if (njsgs(isg) == 3) {
          write(85, format_312), k2sgs(isg, 3), px3, py3, pz3, pmsgs(isg,
            3), inatt, katt(inatt, 1);
        }
      }
      //Clin-5/02/2011
      //C
    }
  }
  //C     number of hadrons formed from partons inside ZPC:
  cmn.nattzp = natt;
  mstj(24) = mstj24;
  //C
}

struct getnp_save
{
  int i;
  int nspec1;
  int nspec2;

  getnp_save() :
    i(fem::int0),
    nspec1(fem::int0),
    nspec2(fem::int0)
  {}
};

//C
//C=======================================================================
void
getnp(
  common& cmn)
{
  FEM_CMN_SVE(getnp);
  // COMMON hmain1
  int& natt = cmn.natt;
  // COMMON hmain2
  const int maxstr = 150001;
  arr_cref<int, 2> katt(cmn.katt, dimension(maxstr, 4));
  arr_cref<float, 2> patt(cmn.patt, dimension(maxstr, 4));
  // COMMON hparnt
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  // COMMON snn
  int& npart1 = cmn.npart1;
  int& npart2 = cmn.npart2;
  float& epsipz = cmn.epsipz;
  float& epsipt = cmn.epsipt;
  float& pzproj = cmn.pzproj;
  float& pztarg = cmn.pztarg;
  //
  // SAVE
  int& i = sve.i;
  int& nspec1 = sve.nspec1;
  int& nspec2 = sve.nspec2;
  //
  //C
  //Cc      SAVE /HMAIN1/
  //Cc      SAVE /HMAIN2/
  //Cc      SAVE /HPARNT/
  //Cc      SAVE /snn/
  //C
  if (natt == 0) {
    npart1 = 0;
    npart2 = 0;
    return;
  }
  //C
  pzproj = fem::sqrt(fem::pow2(hint1(6)) - fem::pow2(hint1(8)));
  pztarg = fem::sqrt(fem::pow2(hint1(7)) - fem::pow2(hint1(9)));
  epsipz = 0.01f;
  //Clin-9/2011-add Pt tolerance in determining spectator nucleons
  //C     (affect string melting runs when LAB frame is used):
  epsipt = 1e-6f;
  //C
  nspec1 = 0;
  nspec2 = 0;
  FEM_DO_SAFE(i, 1, natt) {
    //Clin-9/2011 determine spectator nucleons consistently
    //C           if((KATT(I,1).eq.2112.or.KATT(I,1).eq.2212)
    //C     1          .and.PATT(I, 1).eq.0.and.PATT(I, 2).eq.0) then
    if ((katt(i, 1) == 2112 || katt(i, 1) == 2212) && fem::abs(patt(i,
        1)) <= epsipt && fem::abs(patt(i, 2)) <= epsipt) {
      if (patt(i, 3) > fem::amax1(0.f, pzproj - epsipz)) {
        nspec1++;
      }
      else if (patt(i, 3) < (-pztarg + epsipz)) {
        nspec2++;
      }
    }
  }
  npart1 = ihnt2(1) - nspec1;
  npart2 = ihnt2(3) - nspec2;
  //C
}

struct resdec_save
{
  float dpdecp;
  float enet;
  float esave;
  int idau;
  int ip;
  int irun;
  int kdaut;
  int kf;
  int ksave;
  int lbdaut;
  int ndaut;
  float pxsave;
  float pysave;
  float pzsave;
  float tau0;
  float taudcy;
  float xmsave;

  resdec_save() :
    dpdecp(fem::float0),
    enet(fem::float0),
    esave(fem::float0),
    idau(fem::int0),
    ip(fem::int0),
    irun(fem::int0),
    kdaut(fem::int0),
    kf(fem::int0),
    ksave(fem::int0),
    lbdaut(fem::int0),
    ndaut(fem::int0),
    pxsave(fem::float0),
    pysave(fem::float0),
    pzsave(fem::float0),
    tau0(fem::float0),
    taudcy(fem::float0),
    xmsave(fem::float0)
  {}
};

//C
//C=======================================================================
//C     2/18/03 use PYTHIA to decay eta,rho,omega,k*,phi and Delta
//C     4/2012 added pi0 decay flag:
//C       ipion=0: resonance or pi0 in lb(i1); >0: pi0 in lpion(ipion).
void
resdec(
  common& cmn,
  int const& i1,
  int const& nt,
  int& nnn,
  float const& wid,
  int const& idecay,
  int const& ipion)
{
  FEM_CMN_SVE(resdec);
  common_write write(cmn);
  int& ntmax = cmn.ntmax;
  int& n = cmn.n;
  arr_ref<int, 2> k(cmn.k, dimension(9000, 5));
  arr_ref<float, 2> p(static_cast<common_lujets&>(cmn).p, dimension(9000, 5));
  const int maxstr = 150001;
  arr_ref<float> e(cmn.e, dimension(maxstr));
  arr_ref<int> lb(cmn.lb, dimension(maxstr));
  const int maxr = 1;
  arr_ref<float, 3> rpion(cmn.rpion, dimension(3, maxstr, maxr));
  arr_ref<float, 3> ppion(cmn.ppion, dimension(3, maxstr, maxr));
  arr_ref<float, 2> epion(cmn.epion, dimension(maxstr, maxr));
  arr_ref<int, 2> lpion(cmn.lpion, dimension(maxstr, maxr));
  int& nsav = cmn.nsav;
  int& lb1 = cmn.lb1;
  float& px1 = static_cast<common_leadng&>(cmn).px1;
  float& py1 = static_cast<common_leadng&>(cmn).py1;
  float& pz1 = static_cast<common_leadng&>(cmn).pz1;
  float& em1 = cmn.em1;
  float& e1 = cmn.e1;
  float& xfnl = cmn.xfnl;
  float& yfnl = cmn.yfnl;
  float& zfnl = cmn.zfnl;
  float& tfnl = cmn.tfnl;
  arr_ref<float, 2> tfdpi(cmn.tfdpi, dimension(maxstr, maxr));
  int& nseed = cmn.nseed;
  arr_cref<float> dpertp(cmn.dpertp, dimension(maxstr));
  arr_ref<float, 2> dppion(cmn.dppion, dimension(maxstr, maxr));
  int& ipi0dcy = cmn.ipi0dcy;
  //
  float& dpdecp = sve.dpdecp;
  float& enet = sve.enet;
  float& esave = sve.esave;
  int& idau = sve.idau;
  int& ip = sve.ip;
  int& irun = sve.irun;
  int& kdaut = sve.kdaut;
  int& kf = sve.kf;
  int& ksave = sve.ksave;
  int& lbdaut = sve.lbdaut;
  int& ndaut = sve.ndaut;
  float& pxsave = sve.pxsave;
  float& pysave = sve.pysave;
  float& pzsave = sve.pzsave;
  float& tau0 = sve.tau0;
  float& taudcy = sve.taudcy;
  float& xmsave = sve.xmsave;
  const float apich = 0.140f;
  const float api0 = 0.135f;
  const float addm = 0.02f;
  const float ak0 = 0.498f;
  const float an = 0.940f;
  const float hbarc = 0.19733f;
  //C
  //Cc      SAVE /INPUT2/
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /LUDAT1/
  //Cc      SAVE /LUDAT2/
  //Cc      SAVE /LUDAT3/
  //Cc      SAVE /CC/
  //Cc      SAVE /EE/
  //Cc      SAVE /PA/
  //Cc      SAVE /PB/
  //Cc      SAVE /PC/
  //Cc      SAVE /PD/
  //Cc      SAVE /input1/
  //Cc      SAVE /resdcy/
  //Cc      SAVE /leadng/
  //Cc      SAVE /tdecay/
  //Cc      SAVE /RNDF77/
  irun = idecay;
  //Clin-4/2012 for option of pi0 decay:
  if (nt == ntmax && ipi0dcy == 1 && ((lb1 == 4 && ipion == 0) || ipion >= 1)) {
    kf = 111;
    //C        if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
  }
  else if (lb1 == 0 || lb1 == 25 || lb1 == 26 || lb1 == 27 ||
    lb1 == 28 || lb1 == 29 || fem::iabs(lb1) == 30 || lb1 == 24 || (
    fem::iabs(lb1) >= 6 && fem::iabs(lb1) <= 9) || fem::iabs(
    lb1) == 16) {
    kf = invflv(lb1);
  }
  else {
    return;
  }
  //C
  ip = 1;
  //C     label as undecayed and the only particle in the record:
  n = 1;
  k(ip, 1) = 1;
  k(ip, 3) = 0;
  k(ip, 4) = 0;
  k(ip, 5) = 0;
  //C
  k(ip, 2) = kf;
  //Clin-4/2012 for option of pi0 decay:
  if (ipion == 0) {
    //C
    p(ip, 1) = px1;
    p(ip, 2) = py1;
    p(ip, 3) = pz1;
    //C        em1a=em1
    //C     eta or omega in ART may be below or too close to (pi+pi-pi0) mass,
    //C     causing LUDECY error,thus increase their mass ADDM above this thresh,
    //C     noting that rho (m=0.281) too close to 2pi thrshold fails to decay:
    if ((lb1 == 0 || lb1 == 28) && em1 < (2 * apich + api0 + addm)) {
      em1 = 2 * apich + api0 + addm;
      //C     rho
    }
    else if (lb1 >= 25 && lb1 <= 27 && em1 < (2 * apich + addm)) {
      em1 = 2 * apich + addm;
      //C     K*
    }
    else if (fem::iabs(lb1) == 30 && em1 < (apich + ak0 + addm)) {
      em1 = apich + ak0 + addm;
      //C     Delta created in ART may be below (n+pich) mass, causing LUDECY error:
    }
    else if (fem::iabs(lb1) >= 6 && fem::iabs(lb1) <= 9 && em1 < (
      apich + an + addm)) {
      em1 = apich + an + addm;
    }
    //C        if(em1.ge.(em1a+0.01)) write (6,*)
    //C     1       'Mass increase in resdec():',nt,em1-em1a,lb1
    e1 = fem::sqrt(fem::pow2(em1) + fem::pow2(px1) + fem::pow2(py1) +
      fem::pow2(pz1));
    p(ip, 4) = e1;
    p(ip, 5) = em1;
    //Clin-5/2008:
    dpdecp = dpertp(i1);
    //Clin-4/2012 for option of pi0 decay:
  }
  else if (nt == ntmax && ipi0dcy == 1 && ipion >= 1) {
    p(ip, 1) = ppion(1, ipion, irun);
    p(ip, 2) = ppion(2, ipion, irun);
    p(ip, 3) = ppion(3, ipion, irun);
    p(ip, 5) = epion(ipion, irun);
    p(ip, 4) = fem::sqrt(fem::pow2(p(ip, 5)) + fem::pow2(p(ip, 1)) +
      fem::pow2(p(ip, 2)) + fem::pow2(p(ip, 3)));
    dpdecp = dppion(ipion, irun);
    //Ctest off
    //C           write(99,*) P(IP,4), P(IP,5), dpdecp, ipion, wid
  }
  else {
    write(6, star), "stopped in resdec() a";
    FEM_STOP(0);
  }
  //C
  ludecy(ip);
  //C     add decay time to daughter's formation time at the last timestep:
  if (nt == ntmax) {
    tau0 = hbarc / wid;
    taudcy = tau0 * (-1.f) * fem::alog(1.f - ranart(nseed));
    ndaut = n - nsav;
    if (ndaut <= 1) {
      write(10, star), "note: ndaut(<1)=", ndaut;
      lulist(2);
      FEM_STOP(0);
    }
    //C     lorentz boost:
    //Clin-4/2012 for option of pi0 decay:
    if (ipion == 0) {
      taudcy = taudcy * e1 / em1;
      tfnl += taudcy;
      xfnl += px1 / e1 * taudcy;
      yfnl += py1 / e1 * taudcy;
      zfnl += pz1 / e1 * taudcy;
    }
    else if (ipion >= 1) {
      taudcy = taudcy * p(ip, 4) / p(ip, 5);
      tfnl = tfdpi(ipion, irun) + taudcy;
      xfnl = rpion(1, ipion, irun) + p(ip, 1) / p(ip, 4) * taudcy;
      yfnl = rpion(2, ipion, irun) + p(ip, 2) / p(ip, 4) * taudcy;
      zfnl = rpion(3, ipion, irun) + p(ip, 3) / p(ip, 4) * taudcy;
    }
    else {
      write(6, star), "stopped in resdec() b", ipion, wid, p(ip, 4);
      FEM_STOP(0);
    }
    //C     at the last timestep, assign rho, K0S or eta (decay daughter)
    //C     to lb(i1) only (not to lpion) in order to decay them again:
    //Clin-4/2012 for option of pi0 decay:
    //C           if(n.ge.(nsav+2)) then
    if (n >= (nsav + 2) && ipion == 0) {
      FEM_DO_SAFE(idau, nsav + 2, n) {
        kdaut = k(idau, 2);
        if (kdaut == 221 || kdaut == 113 || kdaut == 213 ||
            kdaut ==  - 213 || kdaut == 310) {
          //C     switch idau and i1(nsav+1):
          ksave = kdaut;
          pxsave = p(idau, 1);
          pysave = p(idau, 2);
          pzsave = p(idau, 3);
          esave = p(idau, 4);
          xmsave = p(idau, 5);
          k(idau, 2) = k(nsav + 1, 2);
          p(idau, 1) = p(nsav + 1, 1);
          p(idau, 2) = p(nsav + 1, 2);
          p(idau, 3) = p(nsav + 1, 3);
          p(idau, 4) = p(nsav + 1, 4);
          p(idau, 5) = p(nsav + 1, 5);
          k(nsav + 1, 2) = ksave;
          p(nsav + 1, 1) = pxsave;
          p(nsav + 1, 2) = pysave;
          p(nsav + 1, 3) = pzsave;
          p(nsav + 1, 4) = esave;
          p(nsav + 1, 5) = xmsave;
          //C     note: phi decay may produce rho, K0s or eta, N*(1535) decay may produce
          //C     eta, but only one daughter may be rho, K0s or eta:
          goto statement_111;
        }
      }
    }
    statement_111:
    //C
    enet = 0.f;
    FEM_DO_SAFE(idau, nsav + 1, n) {
      enet += p(idau, 4);
    }
    //C           if(abs(enet-e1).gt.0.02)
    //C     1          write(93,*) 'resdec(): nt=',nt,enet-e1,lb1
  }
  //C
  FEM_DO_SAFE(idau, nsav + 1, n) {
    kdaut = k(idau, 2);
    lbdaut = iarflv(kdaut);
    //C     K0S and K0L are named K+/K- during hadron cascade, and only
    //C     at the last timestep they keep their real LB # before output;
    //C     K0/K0bar (from K* decay) converted to K0S and K0L at the last timestep:
    if (nt == ntmax && (kdaut == 130 || kdaut == 310 || fem::iabs(
        kdaut) == 311)) {
      if (kdaut == 130) {
        lbdaut = 22;
      }
      else if (kdaut == 310) {
        lbdaut = 24;
      }
      else if (fem::iabs(kdaut) == 311) {
        if (ranart(nseed) < 0.5f) {
          lbdaut = 22;
        }
        else {
          lbdaut = 24;
        }
      }
    }
    //C
    if (idau == (nsav + 1)) {
      //Clin-4/2012 for option of pi0 decay:
      if (ipion == 0) {
        lb(i1) = lbdaut;
        e(i1) = p(idau, 5);
        cmn.px1n = p(idau, 1);
        cmn.py1n = p(idau, 2);
        cmn.pz1n = p(idau, 3);
        //Clin-5/2008:
        cmn.dp1n = dpdecp;
      }
      else if (ipion >= 1) {
        lpion(ipion, irun) = lbdaut;
        epion(ipion, irun) = p(idau, 5);
        ppion(1, ipion, irun) = p(idau, 1);
        ppion(2, ipion, irun) = p(idau, 2);
        ppion(3, ipion, irun) = p(idau, 3);
        rpion(1, ipion, irun) = xfnl;
        rpion(2, ipion, irun) = yfnl;
        rpion(3, ipion, irun) = zfnl;
        tfdpi(ipion, irun) = tfnl;
        dppion(ipion, irun) = dpdecp;
      }
      //C
    }
    else {
      nnn++;
      lpion(nnn, irun) = lbdaut;
      epion(nnn, irun) = p(idau, 5);
      ppion(1, nnn, irun) = p(idau, 1);
      ppion(2, nnn, irun) = p(idau, 2);
      ppion(3, nnn, irun) = p(idau, 3);
      rpion(1, nnn, irun) = xfnl;
      rpion(2, nnn, irun) = yfnl;
      rpion(3, nnn, irun) = zfnl;
      tfdpi(nnn, irun) = tfnl;
      //Clin-5/2008:
      dppion(nnn, irun) = dpdecp;
    }
  }
}

//C
//C=======================================================================
void
inidcy(
  common& cmn)
{
  // COMMON lujets
  int& n = cmn.n;
  //
  //C
  //Cc      SAVE /LUJETS/
  //Cc      SAVE /resdcy/
  n = 1;
  cmn.nsav = n;
}

struct local_save
{
  double detdy;
  double drt;
  double eta0;
  double etcrit;
  double ettest;
  int ip;
  int it;
  int itest;
  double rap0;
  double x0;
  double xtest;
  double y0;
  double ytest;

  local_save() :
    detdy(fem::double0),
    drt(fem::double0),
    eta0(fem::double0),
    etcrit(fem::double0),
    ettest(fem::double0),
    ip(fem::int0),
    it(fem::int0),
    itest(fem::int0),
    rap0(fem::double0),
    x0(fem::double0),
    xtest(fem::double0),
    y0(fem::double0),
    ytest(fem::double0)
  {}
};

//C
//C=======================================================================
//Clin-6/06/02 local parton freezeout motivated from critical density:
void
local(
  common& cmn,
  double const& t)
{
  FEM_CMN_SVE(local);
  common_write write(cmn);
  int& mul = cmn.mul;
  const int maxptn = 400001;
  arr_cref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_cref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_cref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_cref<double> px5(cmn.px5, dimension(maxptn));
  arr_cref<double> py5(cmn.py5, dimension(maxptn));
  arr_cref<double> pz5(cmn.pz5, dimension(maxptn));
  arr_cref<double> e5(cmn.e5, dimension(maxptn));
  arr_cref<double> xmass5(cmn.xmass5, dimension(maxptn));
  arr_cref<int> ityp5(cmn.ityp5, dimension(maxptn));
  arr_ref<double> gxfrz(cmn.gxfrz, dimension(maxptn));
  arr_ref<double> gyfrz(cmn.gyfrz, dimension(maxptn));
  arr_ref<double> gzfrz(cmn.gzfrz, dimension(maxptn));
  arr_ref<double> ftfrz(cmn.ftfrz, dimension(maxptn));
  arr_ref<double> pxfrz(cmn.pxfrz, dimension(maxptn));
  arr_ref<double> pyfrz(cmn.pyfrz, dimension(maxptn));
  arr_ref<double> pzfrz(cmn.pzfrz, dimension(maxptn));
  arr_ref<double> efrz(cmn.efrz, dimension(maxptn));
  arr_ref<double> xmfrz(cmn.xmfrz, dimension(maxptn));
  arr_cref<double> tfrz(cmn.tfrz, dimension(302));
  arr_ref<int> ifrz(cmn.ifrz, dimension(maxptn));
  arr_ref<int> idfrz(cmn.idfrz, dimension(maxptn));
  int& itlast = cmn.itlast;
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  arr_cref<double> eta(cmn.eta, dimension(maxptn));
  arr_cref<double> rap(cmn.rap, dimension(maxptn));
  //
  double& detdy = sve.detdy;
  double& drt = sve.drt;
  double& eta0 = sve.eta0;
  double& etcrit = sve.etcrit;
  double& ettest = sve.ettest;
  int& ip = sve.ip;
  int& it = sve.it;
  int& itest = sve.itest;
  double& rap0 = sve.rap0;
  double& x0 = sve.x0;
  double& xtest = sve.xtest;
  double& y0 = sve.y0;
  double& ytest = sve.ytest;
  const double r0 = 1e0;
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /prec2/
  //Cc      SAVE /frzprc/
  //Cc      SAVE /prec4/
  //Cc      SAVE /prec5/
  //Cc      SAVE /coal/
  //C
  FEM_DO_SAFE(it, 1, 301) {
    if (t >= tfrz(it) && t < tfrz(it + 1)) {
      if (it == itlast) {
        return;
      }
      else {
        itlast = it;
        goto statement_50;
      }
    }
  }
  write(1, star), "local time out of range in LOCAL, stop", t, it;
  FEM_STOP(0);
  statement_50:
  //C
  FEM_DO_SAFE(ip, 1, mul) {
    //C     skip partons which have frozen out:
    if (ifrz(ip) == 1) {
      goto statement_200;
    }
    if (it == 301) {
      //C     freezeout all the left partons beyond the time of 3000 fm:
      etcrit = 1e6;
      goto statement_150;
    }
    else {
      //C     freezeout when transverse energy density < etcrit:
      etcrit = (cmn.ecritl * 2e0 / 3e0);
    }
    //C     skip partons which have not yet formed:
    if (t < ft5(ip)) {
      goto statement_200;
    }
    rap0 = rap(ip);
    eta0 = eta(ip);
    x0 = gx5(ip) + vx(ip) * (t - ft5(ip));
    y0 = gy5(ip) + vy(ip) * (t - ft5(ip));
    detdy = 0e0;
    FEM_DO_SAFE(itest, 1, mul) {
      //C     skip self and partons which have not yet formed:
      if (itest == ip || t < ft5(itest)) {
        goto statement_100;
      }
      ettest = eta(itest);
      xtest = gx5(itest) + vx(itest) * (t - ft5(itest));
      ytest = gy5(itest) + vy(itest) * (t - ft5(itest));
      drt = fem::sqrt(fem::pow2((xtest - x0)) + fem::pow2((ytest - y0)));
      //C     count partons within drt<1 and -1<(eta-eta0)<1:
      if (fem::dabs(ettest - eta0) <= 1e0 && drt <= r0) {
        detdy += fem::dsqrt(fem::pow2(px5(itest)) + fem::pow2(py5(
          itest)) + fem::pow2(xmass5(itest))) * 0.5e0;
      }
      statement_100:;
    }
    detdy = detdy * (fem::pow2(dcosh(eta0))) / (t * 3.1416e0 *
      fem::pow2(r0) * dcosh(rap0));
    //C     when density is below critical density for phase transition, freeze out:
    statement_150:
    if (detdy <= etcrit) {
      ifrz(ip) = 1;
      idfrz(ip) = ityp5(ip);
      pxfrz(ip) = px5(ip);
      pyfrz(ip) = py5(ip);
      pzfrz(ip) = pz5(ip);
      efrz(ip) = e5(ip);
      xmfrz(ip) = xmass5(ip);
      if (t > ft5(ip)) {
        gxfrz(ip) = x0;
        gyfrz(ip) = y0;
        gzfrz(ip) = gz5(ip) + vz(ip) * (t - ft5(ip));
        ftfrz(ip) = t;
      }
      else {
        //C     if this freezeout time < formation time, use formation time & positions.
        //C     This ensures the recovery of default hadron when e_crit=infty:
        gxfrz(ip) = gx5(ip);
        gyfrz(ip) = gy5(ip);
        gzfrz(ip) = gz5(ip);
        ftfrz(ip) = ft5(ip);
      }
    }
    statement_200:;
  }
  //C
}

struct inifrz_save
{
  int it;
  double step1;
  double step2;
  double step3;
  double step4;

  inifrz_save() :
    it(fem::int0),
    step1(fem::double0),
    step2(fem::double0),
    step3(fem::double0),
    step4(fem::double0)
  {}
};

//C
//C=======================================================================
//Clin-6/06/02 initialization for local parton freezeout
void
inifrz(
  common& cmn)
{
  FEM_CMN_SVE(inifrz);
  // COMMON frzprc
  arr_ref<double> tfrz(cmn.tfrz, dimension(302));
  //
  // SAVE
  int& it = sve.it;
  double& step1 = sve.step1;
  double& step2 = sve.step2;
  double& step3 = sve.step3;
  double& step4 = sve.step4;
  //
  //C
  //Cc      SAVE /ilist5/
  //Cc      SAVE /frzprc/
  //C
  //C     for freezeout time 0-10fm, use interval of 0.1fm;
  //C     for 10-100fm, use interval of 1fm;
  //C     for 100-1000fm, use interval of 10fm;
  //C     for 1000-3000fm, use interval of 100fm:
  step1 = 0.1e0;
  step2 = 1e0;
  step3 = 10e0;
  step4 = 100e0;
  //C
  FEM_DO_SAFE(it, 1, 101) {
    tfrz(it) = 0e0 + fem::dble(it - 1) * step1;
  }
  FEM_DO_SAFE(it, 102, 191) {
    tfrz(it) = 10e0 + fem::dble(it - 101) * step2;
  }
  FEM_DO_SAFE(it, 192, 281) {
    tfrz(it) = 100e0 + fem::dble(it - 191) * step3;
  }
  FEM_DO_SAFE(it, 282, 301) {
    tfrz(it) = 1000e0 + fem::dble(it - 281) * step4;
  }
  tfrz(302) = cmn.tlarge;
  //C
}

struct flowh_save
{
  float ene;
  int i;
  int ia;
  int ianh;
  int ic;
  int ifanim;
  int ii;
  int iloop;
  int iy;
  int j;
  int mult;
  int nhadrn;
  double pt2;
  float px;
  float py;
  float rap;
  arr<float> tsh;
  double v2hadr;
  arr<double> v2hevt;
  float xperp2;

  flowh_save() :
    ene(fem::float0),
    i(fem::int0),
    ia(fem::int0),
    ianh(fem::int0),
    ic(fem::int0),
    ifanim(fem::int0),
    ii(fem::int0),
    iloop(fem::int0),
    iy(fem::int0),
    j(fem::int0),
    mult(fem::int0),
    nhadrn(fem::int0),
    pt2(fem::double0),
    px(fem::float0),
    py(fem::float0),
    rap(fem::float0),
    tsh(dimension(31), fem::fill0),
    v2hadr(fem::double0),
    v2hevt(dimension(3), fem::fill0),
    xperp2(fem::float0)
  {}
};

//C
//C$$$clin-5/2009 v2 analysis
//C$$$c=======================================================================
//C$$$c     idd=0,1,2,3 specifies different subroutines for partonic flow analysis.
//C$$$        subroutine flowp(idd)
//C$$$c
//C$$$        implicit double precision  (a-h, o-z)
//C$$$        real dt
//C$$$        parameter (MAXPTN=400001)
//C$$$csp
//C$$$        parameter (bmt=0.05d0)
//C$$$        dimension nlfile(3),nsfile(3),nmfile(3)
//C$$$c
//C$$$        dimension v2pp(3),xnpp(3),v2psum(3),v2p2sm(3),nfile(3)
//C$$$        dimension tsp(31),v2pevt(3),v2pavg(3),varv2p(3)
//C$$$        common /ilist1/
//C$$$     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
//C$$$     &     ictype, icsta(MAXPTN),
//C$$$     &     nic(MAXPTN), icels(MAXPTN)
//C$$$cc      SAVE /ilist1/
//C$$$        COMMON /para1/ mul
//C$$$cc      SAVE /para1/
//C$$$        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
//C$$$     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
//C$$$     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
//C$$$cc      SAVE /prec2/
//C$$$        COMMON /pflow/ v2p(30,3),xnpart(30,3),etp(30,3),
//C$$$     1       s2p(30,3),v2p2(30,3),nevt(30)
//C$$$cc      SAVE /pflow/
//C$$$        COMMON /pflowf/ v2pf(30,3),xnpf(30,3),etpf(30,3),
//C$$$     1                 xncoll(30),s2pf(30,3),v2pf2(30,3)
//C$$$cc      SAVE /pflowf/
//C$$$        COMMON /pfrz/ v2pfrz(30,3),xnpfrz(30,3),etpfrz(30,3),
//C$$$     1       s2pfrz(30,3),v2p2fz(30,3),tscatt(31),
//C$$$     2       nevtfz(30),iscatt(30)
//C$$$cc      SAVE /pfrz/
//C$$$        COMMON /hflow/ v2h(30,3),xnhadr(30,3),eth(30,3),
//C$$$     1 v2h2(30,3),s2h(30,3)
//C$$$cc      SAVE /hflow/
//C$$$        COMMON /AREVT/ IAEVT, IARUN, MISS
//C$$$cc      SAVE /AREVT/
//C$$$        common/anim/nevent,isoft,isflag,izpc
//C$$$cc      SAVE /anim/
//C$$$        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
//C$$$cc      SAVE /input1/
//C$$$        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE,
//C$$$     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
//C$$$cc      SAVE /INPUT2/
//C$$$cc      SAVE itimep,iaevtp,v2pp,xnpp,v2psum,v2p2sm
//C$$$cc      SAVE nfile,itanim,nlfile,nsfile,nmfile
//C$$$        common /precpb/vxp(MAXPTN),vyp(MAXPTN),vzp(MAXPTN)
//C$$$        SAVE
//C$$$csp
//C$$$        dimension etpl(30,3),etps(30,3),etplf(30,3),etpsf(30,3),
//C$$$     &       etlfrz(30,3),etsfrz(30,3),
//C$$$     &       xnpl(30,3),xnps(30,3),xnplf(30,3),xnpsf(30,3),
//C$$$     &       xnlfrz(30,3),xnsfrz(30,3),
//C$$$     &       v2pl(30,3),v2ps(30,3),v2plf(30,3),v2psf(30,3),
//C$$$     &       s2pl(30,3),s2ps(30,3),s2plf(30,3),s2psf(30,3),
//C$$$     &       DMYil(50,3),DMYfl(50,3),
//C$$$     &       DMYis(50,3),DMYfs(50,3)
//C$$$        data tsp/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
//C$$$     &       1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
//C$$$     &       2  , 3,   4,   5,   6,   7,   8,   9,   10,  20,  30/
//C$$$c     idd=0: initialization for flow analysis, called by artdri.f:
//C$$$        if(idd.eq.0) then
//C$$$           nfile(1)=60
//C$$$           nfile(2)=64
//C$$$           nfile(3)=20
//C$$$           OPEN (nfile(1),FILE='ana1/v2p.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(1)+1,
//C$$$     1 FILE = 'ana1/v2p-formed.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(1)+2,
//C$$$     1 FILE = 'ana1/v2p-active.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(1)+3,
//C$$$     1 FILE = 'ana1/v2ph.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(2),FILE='ana1/v2p-y2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(2)+1,
//C$$$     1 FILE = 'ana1/v2p-formed2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(2)+2,
//C$$$     1 FILE = 'ana1/v2p-active2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(2)+3,
//C$$$     1 FILE = 'ana1/v2ph-y2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(3),FILE='ana1/v2p-y1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(3)+1,
//C$$$     1 FILE = 'ana1/v2p-formed1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(3)+2,
//C$$$     1 FILE = 'ana1/v2p-active1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nfile(3)+3,
//C$$$     1 FILE = 'ana1/v2ph-y1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (49, FILE = 'ana1/v2p-ebe.dat', STATUS = 'UNKNOWN')
//C$$$           write(49, *) '    ievt,  v2p,  v2p_y2,   v2p_y1'
//C$$$c
//C$$$           OPEN (59, FILE = 'ana1/v2h.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (68, FILE = 'ana1/v2h-y2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (69, FILE = 'ana1/v2h-y1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (88, FILE = 'ana1/v2h-ebe.dat', STATUS = 'UNKNOWN')
//C$$$           write(88, *) '    ievt,  v2h,  v2h_y2,   v2h_y1'
//C$$$csp07/05
//C$$$           nlfile(1)=70
//C$$$           nlfile(2)=72
//C$$$           nlfile(3)=74
//C$$$           OPEN (nlfile(1),FILE='ana1/mtl.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nlfile(1)+1,
//C$$$     1 FILE = 'ana1/mtl-formed.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nlfile(2),FILE='ana1/mtl-y2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nlfile(2)+1,
//C$$$     1 FILE = 'ana1/mtl-formed2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nlfile(3),FILE='ana1/mtl-y1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nlfile(3)+1,
//C$$$     1 FILE = 'ana1/mtl-formed1.dat', STATUS = 'UNKNOWN')
//C$$$           nsfile(1)=76
//C$$$           nsfile(2)=78
//C$$$           nsfile(3)=80
//C$$$           OPEN (nsfile(1),FILE='ana1/mts.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nsfile(1)+1,
//C$$$     1 FILE = 'ana1/mts-formed.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nsfile(2),FILE='ana1/mts-y2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nsfile(2)+1,
//C$$$     1 FILE = 'ana1/mts-formed2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nsfile(3),FILE='ana1/mts-y1.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nsfile(3)+1,
//C$$$     1 FILE = 'ana1/mts-formed1.dat', STATUS = 'UNKNOWN')
//C$$$           nmfile(1)=82
//C$$$           nmfile(2)=83
//C$$$           nmfile(3)=84
//C$$$           OPEN (nmfile(1),FILE='ana1/Nmt.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nmfile(2),FILE='ana1/Nmt-y2.dat', STATUS = 'UNKNOWN')
//C$$$           OPEN (nmfile(3),FILE='ana1/Nmt-y1.dat', STATUS = 'UNKNOWN')
//C$$$clin-8/2015: changed unit number of animation files,
//C$$$ctest off     turn off animation output (0 to turn off and 1 to turn on):
//C$$$           ifanim=0
//C$$$clin-11/27/00 for animation:
//C$$$           if(ifanim.eq.1) then
//C$$$              OPEN (10, FILE = 'ana1/h-animate.dat', STATUS = 'UNKNOWN')
//C$$$              write(10,*) ntmax, dt
//C$$$              OPEN (11, FILE = 'ana1/p-animate.dat', STATUS = 'UNKNOWN')
//C$$$              OPEN (15, FILE = 'ana1/p-finalft.dat', STATUS = 'UNKNOWN')
//C$$$           endif
//C$$$clin-10/2014: write out partons at all eta, turn off now:
//C$$$c           if(nevent.ge.1)
//C$$$           if(nevent.lt.1)
//C$$$     1          OPEN (93, FILE = 'ana1/parton-t.dat', STATUS='UNKNOWN')
//C$$$c
//C$$$           itimep=0
//C$$$           itanim=0
//C$$$           iaevtp=0
//C$$$csp
//C$$$           do 1002 ii=1,50
//C$$$              do 1001 iy=1,3
//C$$$                 DMYil(ii,iy) = 0d0
//C$$$                 DMYfl(ii,iy) = 0d0
//C$$$                 DMYis(ii,iy) = 0d0
//C$$$                 DMYfs(ii,iy) = 0d0
//C$$$ 1001         continue
//C$$$ 1002      continue
//C$$$c
//C$$$           do 1003 ii=1,31
//C$$$              tscatt(ii)=0d0
//C$$$ 1003      continue
//C$$$           do 1005 ii=1,30
//C$$$              nevt(ii)=0
//C$$$              xncoll(ii)=0d0
//C$$$              nevtfz(ii)=0
//C$$$              iscatt(ii)=0
//C$$$              do 1004 iy=1,3
//C$$$                 v2p(ii,iy)=0d0
//C$$$                 v2p2(ii,iy)=0d0
//C$$$                 s2p(ii,iy)=0d0
//C$$$                 etp(ii,iy)=0d0
//C$$$                 xnpart(ii,iy)=0d0
//C$$$                 v2pf(ii,iy)=0d0
//C$$$                 v2pf2(ii,iy)=0d0
//C$$$                 s2pf(ii,iy)=0d0
//C$$$                 etpf(ii,iy)=0d0
//C$$$                 xnpf(ii,iy)=0d0
//C$$$                 v2pfrz(ii,iy)=0d0
//C$$$                 v2p2fz(ii,iy)=0d0
//C$$$                 s2pfrz(ii,iy)=0d0
//C$$$                 etpfrz(ii,iy)=0d0
//C$$$                 xnpfrz(ii,iy)=0d0
//C$$$csp07/05
//C$$$                 etpl(ii,iy)=0d0
//C$$$                 etps(ii,iy)=0d0
//C$$$                 etplf(ii,iy)=0d0
//C$$$                 etpsf(ii,iy)=0d0
//C$$$                 etlfrz(ii,iy)=0d0
//C$$$                 etsfrz(ii,iy)=0d0
//C$$$              xnpl(ii,iy)=0d0
//C$$$              xnps(ii,iy)=0d0
//C$$$              xnplf(ii,iy)=0d0
//C$$$              xnpsf(ii,iy)=0d0
//C$$$              xnlfrz(ii,iy)=0d0
//C$$$              xnsfrz(ii,iy)=0d0
//C$$$              v2pl(ii,iy)=0d0
//C$$$              v2ps(ii,iy)=0d0
//C$$$              v2plf(ii,iy)=0d0
//C$$$              v2psf(ii,iy)=0d0
//C$$$              s2pl(ii,iy)=0d0
//C$$$              s2ps(ii,iy)=0d0
//C$$$              s2plf(ii,iy)=0d0
//C$$$              s2psf(ii,iy)=0d0
//C$$$ 1004      continue
//C$$$ 1005   continue
//C$$$           do 1006 iy=1,3
//C$$$              v2pevt(iy)=0d0
//C$$$              v2pavg(iy)=0d0
//C$$$              varv2p(iy)=0d0
//C$$$              v2pp(iy)=0.d0
//C$$$              xnpp(iy)=0d0
//C$$$              v2psum(iy)=0.d0
//C$$$              v2p2sm(iy)=0.d0
//C$$$ 1006      continue
//C$$$c     idd=1: calculate parton elliptic flow, called by zpc.f:
//C$$$        else if(idd.eq.1) then
//C$$$           if(iaevt.ne.iaevtp.and.ianp.eq.31) itanim=0
//C$$$c
//C$$$           t2time = FT5(iscat)
//C$$$           do 1008 ianp = 1, 30
//C$$$              if (t2time.lt.tsp(ianp+1).and.t2time.ge.tsp(ianp)) then
//C$$$c     write flow info only once at each fixed time:
//C$$$                 xncoll(ianp)=xncoll(ianp)+1d0
//C$$$c     to prevent an earlier t2time comes later in the same event
//C$$$c     and mess up nevt:
//C$$$                 if(ianp.le.itimep.and.iaevt.eq.iaevtp) goto 101
//C$$$                 nevt(ianp)=nevt(ianp)+1
//C$$$                 tscatt(ianp+1)=t2time
//C$$$                 iscatt(ianp)=1
//C$$$                 nevtfz(ianp)=nevtfz(ianp)+1
//C$$$                 do 100 I=1,mul
//C$$$clin-8/2015 to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID:
//C$$$c                    ypartn=0.5d0*dlog((E5(i)+PZ5(i))
//C$$$c     1                   /(E5(i)-PZ5(i)+1.d-8))
//C$$$                    delta=1d-8
//C$$$                    if((E5(i)-dabs(PZ5(i))+delta).le.0) then
//C$$$                       ypartn=1000000.d0*sign(1.d0,PZ5(i))
//C$$$                       write(6,*) 'ypartn error',E5(i)-dabs(PZ5(i))
//C$$$                    else
//C$$$                       ypartn=0.5d0*dlog((E5(i)+PZ5(i)+delta)
//C$$$     1                      /(E5(i)-PZ5(i)+delta))
//C$$$                    endif
//C$$$                    pt2=PX5(I)**2+PY5(I)**2
//C$$$ctest off: initial (pt,y) and (x,y) distribution:
//C$$$c                    idtime=1
//C$$$c                    if(ianp.eq.idtime) then
//C$$$c                       iityp=iabs(ITYP5(I))
//C$$$c                       if(iityp.eq.1.or.iityp.eq.2) then
//C$$$c                          write(651,*) dsqrt(pt2),ypartn
//C$$$c                          write(654,*) GX5(I),GY5(I)
//C$$$c                       elseif(iityp.eq.1103.or.iityp.eq.2101
//C$$$c     1 .or.iityp.eq.2103.or.iityp.eq.2203.
//C$$$c     2 .or.iityp.eq.3101.or.iityp.eq.3103.
//C$$$c     3 .or.iityp.eq.3201.or.iityp.eq.3203.or.iityp.eq.3303)
//C$$$c     4 then
//C$$$c                          write(652,*) dsqrt(pt2),ypartn
//C$$$c                          write(655,*) GX5(I),GY5(I)
//C$$$c                       elseif(iityp.eq.21) then
//C$$$c                          write(653,*) dsqrt(pt2),ypartn
//C$$$c                          write(656,*) GX5(I),GY5(I)
//C$$$c                       endif
//C$$$c                    endif
//C$$$ctest-end
//C$$$ctest off density with 2fm radius and z:(-0.1*t,0.1*t):
//C$$$c                    gx_now=GX5(i)+(t2time-FT5(i))*PX5(i)/E5(i)
//C$$$c                    gy_now=GY5(i)+(t2time-FT5(i))*PY5(i)/E5(i)
//C$$$c                    gz_now=GZ5(i)+(t2time-FT5(i))*PZ5(i)/E5(i)
//C$$$c                    rt_now=dsqrt(gx_now**2+gy_now**2)
//C$$$c                    zmax=0.1d0*t2time
//C$$$c                    volume=3.1416d0*(2d0**2)*(2*zmax)
//C$$$c                    if(rt_now.gt.2d0.or.dabs(gz_now).gt.zmax)
//C$$$c     1                   goto 100
//C$$$ctest-end
//C$$$                    iloop=1
//C$$$                    if(dabs(ypartn).le.1d0) then
//C$$$                       iloop=2
//C$$$                       if(dabs(ypartn).le.0.5d0) then
//C$$$                          iloop=3
//C$$$                       endif
//C$$$                    endif
//C$$$                    do 50 iy=1,iloop
//C$$$clin-5/2012:
//C$$$c                       if(pt2.gt.0.) then
//C$$$                       if(pt2.gt.0d0) then
//C$$$                          v2prtn=(PX5(I)**2-PY5(I)**2)/pt2
//C$$$clin-5/2012:
//C$$$c                          if(abs(v2prtn).gt.1.)
//C$$$                          if(dabs(v2prtn).gt.1d0)
//C$$$     1 write(nfile(iy),*) 'v2prtn>1',v2prtn
//C$$$                          v2p(ianp,iy)=v2p(ianp,iy)+v2prtn
//C$$$                          v2p2(ianp,iy)=v2p2(ianp,iy)+v2prtn**2
//C$$$                       endif
//C$$$                       xperp2=GX5(I)**2+GY5(I)**2
//C$$$clin-5/2012:
//C$$$c                       if(xperp2.gt.0.)
//C$$$                       if(xperp2.gt.0d0)
//C$$$     1        s2p(ianp,iy)=s2p(ianp,iy)+(GX5(I)**2-GY5(I)**2)/xperp2
//C$$$                       xnpart(ianp,iy)=xnpart(ianp,iy)+1d0
//C$$$                       etp(ianp,iy)=etp(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
//C$$$ctest off density:
//C$$$c                       etp(ianp,iy)=etp(ianp,iy)
//C$$$c     1                  +dsqrt(pt2+XMASS5(I)**2+PZ5(i)**2)/volume
//C$$$clin-2/22/00 to write out parton info only for formed ones:
//C$$$                       if(FT5(I).le.t2time) then
//C$$$                          v2pf(ianp,iy)=v2pf(ianp,iy)+v2prtn
//C$$$                          v2pf2(ianp,iy)=v2pf2(ianp,iy)+v2prtn**2
//C$$$clin-5/2012:
//C$$$c                          if(xperp2.gt.0.)
//C$$$                          if(xperp2.gt.0d0)
//C$$$     1        s2pf(ianp,iy)=s2pf(ianp,iy)+(GX5(I)**2-GY5(I)**2)/xperp2
//C$$$                          xnpf(ianp,iy)=xnpf(ianp,iy)+1d0
//C$$$                  etpf(ianp,iy)=etpf(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
//C$$$ctest off density:
//C$$$c                  etpf(ianp,iy)=etpf(ianp,iy)
//C$$$c     1                   +dsqrt(pt2+XMASS5(I)**2+PZ5(i)**2)/volume
//C$$$                       endif
//C$$$ 50                    continue
//C$$$ 100                 continue
//C$$$                 itimep=ianp
//C$$$                 iaevtp=iaevt
//C$$$clin-3/30/00 ebe v2 variables:
//C$$$                 if(ianp.eq.30) then
//C$$$                    do 1007 iy=1,3
//C$$$                       npartn=IDINT(xnpart(ianp,iy)-xnpp(iy))
//C$$$                       if(npartn.ne.0) then
//C$$$                          v2pevt(iy)=(v2p(ianp,iy)-v2pp(iy))/npartn
//C$$$                          v2psum(iy)=v2psum(iy)+v2pevt(iy)
//C$$$                          v2p2sm(iy)=v2p2sm(iy)+v2pevt(iy)**2
//C$$$                          v2pp(iy)=v2p(ianp,iy)
//C$$$                          xnpp(iy)=xnpart(ianp,iy)
//C$$$                       endif
//C$$$ 1007               continue
//C$$$                    write(49, 160) iaevt,v2pevt
//C$$$                 endif
//C$$$                 goto 101
//C$$$              endif
//C$$$ 1008      continue
//C$$$clin-11/28/00 for animation:
//C$$$clin-8/2015 ctest off turn off parton-t.dat:
//C$$$c 101       if(nevent.ge.1) then
//C$$$ 101       if(nevent.lt.1) then
//C$$$              do 110 nt = 1, ntmax
//C$$$                 time1=dble(nt*dt)
//C$$$                 time2=dble((nt+1)*dt)
//C$$$                 if (t2time.lt.time2.and.t2time.ge.time1) then
//C$$$                    if(nt.le.itanim) return
//C$$$                    if(ifanim.eq.1) write(11,*) t2time
//C$$$                    iform=0
//C$$$                    ne1all=0
//C$$$                    ne1form=0
//C$$$                    do 1009 I=1,mul
//C$$$c     Calculate parton coordinates after propagation to current time:
//C$$$                       gz_now=GZ5(i)+(t2time-FT5(i))*PZ5(i)/E5(i)
//C$$$                       If(dabs(gz_now).lt.t2time) then
//C$$$                       etap=0.5d0*dlog((t2time+gz_now)/(t2time-gz_now))
//C$$$                       else
//C$$$                          etap=1000000.d0*sign(1.d0,gz_now)
//C$$$                       endif
//C$$$                       ne1all=ne1all+1
//C$$$                       if(FT5(I).le.t2time) ne1form=ne1form+1
//C$$$c     write out parton info only for formed ones for animation:
//C$$$                       if(FT5(I).le.t2time) iform=iform+1
//C$$$ 1009               continue
//C$$$clin-8/2015 for animation:
//C$$$                    if(ifanim.eq.1) write(11,*) iform
//C$$$                    write(93,184) 'evt#,t,np,npformed=',
//C$$$     1                   iaevt,t2time,ne1all,ne1form
//C$$$ 184                format(a20,i7,f8.4,2(1x,i6))
//C$$$c
//C$$$                    do 120 I=1,mul
//C$$$                       if(FT5(I).le.t2time) then
//C$$$c     propagate formed partons to current time t2time using parton v:
//C$$$                          gz_now=GZ5(i)+(t2time-FT5(i))*PZ5(i)/E5(i)
//C$$$                       else
//C$$$c     back-propagate unformed partons using parent hadron v:
//C$$$                          gz_now=GZ5(i)+(t2time-FT5(i))*vzp(i)
//C$$$                       endif
//C$$$c
//C$$$                       If(dabs(gz_now).lt.t2time) then
//C$$$                       etap=0.5d0*dlog((t2time+gz_now)/(t2time-gz_now))
//C$$$                       else
//C$$$                          etap=1000000.d0*sign(1.d0,gz_now)
//C$$$                       endif
//C$$$c     calculate other coordinates of the parton:
//C$$$                       if(FT5(I).le.t2time) then
//C$$$                          gx_now=GX5(i)+(t2time-FT5(i))*PX5(i)/E5(i)
//C$$$                          gy_now=GY5(i)+(t2time-FT5(i))*PY5(i)/E5(i)
//C$$$                       else
//C$$$                          gx_now=GX5(i)+(t2time-FT5(i))*vxp(i)
//C$$$                          gy_now=GY5(i)+(t2time-FT5(i))*vyp(i)
//C$$$                       endif
//C$$$                       write(93,185) ITYP5(i),PX5(i),PY5(i),PZ5(i),
//C$$$     1                      XMASS5(i),gx_now,gy_now,ft5(i),etap
//C$$$clin-8/2015 for animation:
//C$$$                       if(ifanim.eq.1.and.FT5(I).le.t2time) then
//C$$$                          write(11,180) ITYP5(i),GX5(i),GY5(i),GZ5(i),
//C$$$     1                         FT5(i),PX5(i),PY5(i),PZ5(i),E5(i)
//C$$$                       endif
//C$$$ 185           format(i3,3(1x,f8.3),1x,f8.4,1x,2(f8.3,1x),f11.4,1x,f8.3)
//C$$$ 120                continue
//C$$$                    itanim=nt
//C$$$                 endif
//C$$$ 110          continue
//C$$$           endif
//C$$$c
//C$$$ 160       format(i10,3(2x,f9.5))
//C$$$ 180       format(i6,8(1x,f7.2))
//C$$$clin-5/17/01 calculate v2 for active partons (which have not frozen out):
//C$$$c     idd=3, called at end of zpc.f:
//C$$$        else if(idd.eq.3) then
//C$$$           do 1010 ianp=1,30
//C$$$              if(iscatt(ianp).eq.0) tscatt(ianp+1)=tscatt(ianp)
//C$$$ 1010      continue
//C$$$           do 350 I=1,mul
//C$$$clin-8/2015 to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID:
//C$$$c              ypartn=0.5d0*dlog((E5(i)+PZ5(i))
//C$$$c     1             /(E5(i)-PZ5(i)+1.d-8))
//C$$$              delta=1d-8
//C$$$              if((E5(i)-dabs(PZ5(i))+delta).le.0) then
//C$$$                 write(6,*) 'ypartn error',E5(i)-dabs(PZ5(i))
//C$$$                 ypartn=1000000.d0*sign(1.d0,PZ5(i))
//C$$$              else
//C$$$                 ypartn=0.5d0*dlog((E5(i)+PZ5(i)+delta)
//C$$$     1                /(E5(i)-PZ5(i)+delta))
//C$$$              endif
//C$$$              pt2=PX5(I)**2+PY5(I)**2
//C$$$              iloop=1
//C$$$              if(dabs(ypartn).le.1d0) then
//C$$$                 iloop=2
//C$$$                 if(dabs(ypartn).le.0.5d0) then
//C$$$                    iloop=3
//C$$$                 endif
//C$$$              endif
//C$$$c
//C$$$              do 325 ianp=1,30
//C$$$                 if(iscatt(ianp).ne.0) then
//C$$$                    if(FT5(I).lt.tscatt(ianp+1)
//C$$$     1 .and.FT5(I).ge.tscatt(ianp)) then
//C$$$                       do 1011 iy=1,iloop
//C$$$clin-5/2012:
//C$$$c                          if(pt2.gt.0.) then
//C$$$                          if(pt2.gt.0d0) then
//C$$$                             v2prtn=(PX5(I)**2-PY5(I)**2)/pt2
//C$$$                             v2pfrz(ianp,iy)=v2pfrz(ianp,iy)+v2prtn
//C$$$                     v2p2fz(ianp,iy)=v2p2fz(ianp,iy)+v2prtn**2
//C$$$                          endif
//C$$$                          xperp2=GX5(I)**2+GY5(I)**2
//C$$$clin-5/2012:
//C$$$c                          if(xperp2.gt.0.) s2pfrz(ianp,iy)=
//C$$$                          if(xperp2.gt.0d0) s2pfrz(ianp,iy)=
//C$$$     1 s2pfrz(ianp,iy)+(GX5(I)**2-GY5(I)**2)/xperp2
//C$$$        etpfrz(ianp,iy)=etpfrz(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
//C$$$                          xnpfrz(ianp,iy)=xnpfrz(ianp,iy)+1d0
//C$$$ctest off density:
//C$$$c                    etpfrz(ianp,iy)=etpfrz(ianp,iy)
//C$$$c     1                   +dsqrt(pt2+XMASS5(I)**2+PZ5(i)**2)/volume
//C$$$csp07/05
//C$$$            if(ITYP5(I).eq.1.or.ITYP5(I).eq.2)then
//C$$$              etlfrz(ianp,iy)=etlfrz(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
//C$$$              xnlfrz(ianp,iy)=xnlfrz(ianp,iy)+1d0
//C$$$            elseif(ITYP5(I).eq.3)then
//C$$$              etsfrz(ianp,iy)=etsfrz(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
//C$$$              xnsfrz(ianp,iy)=xnsfrz(ianp,iy)+1d0
//C$$$            endif
//C$$$csp07/05 end
//C$$$ 1011    continue
//C$$$c     parton freezeout info taken, proceed to next parton:
//C$$$                       goto 350
//C$$$                    endif
//C$$$                 endif
//C$$$ 325          continue
//C$$$ 350       continue
//C$$$c     idd=2: calculate average partonic elliptic flow, called from artdri.f,
//C$$$        else if(idd.eq.2) then
//C$$$           do 1012 i=1,3
//C$$$              write(nfile(i),*) '   tsp,   v2p,     v2p2, '
//C$$$     1 '   s2p,  etp,   xmult,    nevt,  xnctot'
//C$$$              write ((nfile(i)+1),*) '  tsp,   v2pf,   v2pf2, '
//C$$$     1 '   s2pf, etpf,  xnform,  xcrate'
//C$$$              write ((nfile(i)+2),*) '  tsp,   v2pa,   v2pa2, '
//C$$$     1 '   s2pa, etpa,  xmult_ap,  xnform,   nevt'
//C$$$              write ((nfile(i)+3),*) '  tsph,  v2ph,   v2ph2, '
//C$$$     1 '   s2ph, etph,  xmult_(ap/2+h),xmult_ap/2,nevt'
//C$$$csp
//C$$$           write(nlfile(i),*) '   tsp,    v2,     s2,    etp,    xmul'
//C$$$           write(nsfile(i),*) '   tsp,    v2,     s2,    etp,    xmul'
//C$$$           write(nlfile(i)+1,*) '   tsp,    v2,     s2,    etp,    xmul'
//C$$$           write(nsfile(i)+1,*) '   tsp,    v2,     s2,    etp,    xmul'
//C$$$c
//C$$$ 1012   continue
//C$$$clin-3/30/00 ensemble average & variance of v2 (over particles in all events):
//C$$$           do 150 ii=1, 30
//C$$$              if(nevt(ii).eq.0) goto 150
//C$$$              do 1014 iy=1,3
//C$$$clin-5/2012:
//C$$$c                 if(xnpart(ii,iy).gt.1) then
//C$$$                 if(xnpart(ii,iy).gt.1d0) then
//C$$$                    v2p(ii,iy)=v2p(ii,iy)/xnpart(ii,iy)
//C$$$                    v2p2(ii,iy)=dsqrt((v2p2(ii,iy)/xnpart(ii,iy)
//C$$$     1                    -v2p(ii,iy)**2)/(xnpart(ii,iy)-1))
//C$$$                    s2p(ii,iy)=s2p(ii,iy)/xnpart(ii,iy)
//C$$$c xmult and etp are multiplicity and et for an averaged event:
//C$$$                    xmult=dble(xnpart(ii,iy)/dble(nevt(ii)))
//C$$$                    etp(ii,iy)=etp(ii,iy)/dble(nevt(ii))
//C$$$csp
//C$$$                    etpl(ii,iy)=etpl(ii,iy)/dble(nevt(ii))
//C$$$                    etps(ii,iy)=etps(ii,iy)/dble(nevt(ii))
//C$$$c
//C$$$                    xnctot=0d0
//C$$$                    do 1013 inum=1,ii
//C$$$                       xnctot=xnctot+xncoll(inum)
//C$$$ 1013               continue
//C$$$                    if(nevt(1).ne.0) xnctot=xnctot/nevt(1)
//C$$$                    write (nfile(iy),200) tsp(ii),v2p(ii,iy),
//C$$$     1      v2p2(ii,iy),s2p(ii,iy),etp(ii,iy),xmult,nevt(ii),xnctot
//C$$$                 endif
//C$$$                 if(nevt(ii).ne.0)
//C$$$     1                xcrate=xncoll(ii)/(tsp(ii+1)-tsp(ii))/nevt(ii)
//C$$$c
//C$$$clin-5/2012:
//C$$$c                 if(xnpf(ii,iy).gt.1) then
//C$$$                 if(xnpf(ii,iy).gt.1d0) then
//C$$$                    v2pf(ii,iy)=v2pf(ii,iy)/xnpf(ii,iy)
//C$$$                    v2pf2(ii,iy)=dsqrt((v2pf2(ii,iy)/xnpf(ii,iy)
//C$$$     1                    -v2pf(ii,iy)**2)/(xnpf(ii,iy)-1))
//C$$$                    s2pf(ii,iy)=s2pf(ii,iy)/xnpf(ii,iy)
//C$$$                    xnform=dble(xnpf(ii,iy)/dble(nevt(ii)))
//C$$$                    etpf(ii,iy)=etpf(ii,iy)/dble(nevt(ii))
//C$$$csp
//C$$$                    etplf(ii,iy)=etplf(ii,iy)/dble(nevt(ii))
//C$$$                    etpsf(ii,iy)=etpsf(ii,iy)/dble(nevt(ii))
//C$$$c
//C$$$                    write (nfile(iy)+1, 210) tsp(ii),v2pf(ii,iy),
//C$$$     1      v2pf2(ii,iy),s2pf(ii,iy),etpf(ii,iy),xnform,xcrate
//C$$$                 endif
//C$$$csp
//C$$$clin-5/2012:
//C$$$c                 if(xnpl(ii,iy).gt.1) then
//C$$$                 if(xnpl(ii,iy).gt.1d0) then
//C$$$                    v2pl(ii,iy)=v2pl(ii,iy)/xnpl(ii,iy)
//C$$$                    s2pl(ii,iy)=s2pl(ii,iy)/xnpl(ii,iy)
//C$$$                    xmult=dble(xnpl(ii,iy)/dble(nevt(ii)))
//C$$$                    etpl(ii,iy)=etpl(ii,iy)/dble(nevt(ii))
//C$$$                    write (nlfile(iy),201) tsp(ii),v2pl(ii,iy),
//C$$$     1        s2pl(ii,iy),etpl(ii,iy),xmult
//C$$$                 endif
//C$$$clin-5/2012:
//C$$$c                 if(xnps(ii,iy).gt.1) then
//C$$$                 if(xnps(ii,iy).gt.1d0) then
//C$$$                    v2ps(ii,iy)=v2ps(ii,iy)/xnps(ii,iy)
//C$$$                    s2ps(ii,iy)=s2ps(ii,iy)/xnps(ii,iy)
//C$$$                    xmult=dble(xnps(ii,iy)/dble(nevt(ii)))
//C$$$                    etps(ii,iy)=etps(ii,iy)/dble(nevt(ii))
//C$$$                    write (nsfile(iy),201) tsp(ii),v2ps(ii,iy),
//C$$$     1        s2ps(ii,iy),etps(ii,iy),xmult
//C$$$                 endif
//C$$$clin-5/2012:
//C$$$c                 if(xnplf(ii,iy).gt.1) then
//C$$$                 if(xnplf(ii,iy).gt.1d0) then
//C$$$                    v2plf(ii,iy)=v2plf(ii,iy)/xnplf(ii,iy)
//C$$$                    s2plf(ii,iy)=s2plf(ii,iy)/xnplf(ii,iy)
//C$$$                    xmult=dble(xnplf(ii,iy)/dble(nevt(ii)))
//C$$$                    etplf(ii,iy)=etplf(ii,iy)/dble(nevt(ii))
//C$$$                    write (nlfile(iy)+1,201) tsp(ii),v2plf(ii,iy),
//C$$$     1        s2plf(ii,iy),etplf(ii,iy),xmult
//C$$$                 endif
//C$$$clin-5/2012:
//C$$$c                 if(xnpsf(ii,iy).gt.1) then
//C$$$                 if(xnpsf(ii,iy).gt.1d0) then
//C$$$                    v2psf(ii,iy)=v2psf(ii,iy)/xnpsf(ii,iy)
//C$$$                    s2psf(ii,iy)=s2psf(ii,iy)/xnpsf(ii,iy)
//C$$$                    xmult=dble(xnpsf(ii,iy)/dble(nevt(ii)))
//C$$$                    etpsf(ii,iy)=etpsf(ii,iy)/dble(nevt(ii))
//C$$$                    write (nsfile(iy)+1,201) tsp(ii),v2psf(ii,iy),
//C$$$     1        s2psf(ii,iy),etpsf(ii,iy),xmult
//C$$$                 endif
//C$$$csp-end
//C$$$ 1014         continue
//C$$$ 150           continue
//C$$$csp07/05 initial & final mt distrb
//C$$$               scalei=0d0
//C$$$               scalef=0d0
//C$$$               if(nevt(1).ne.0) SCALEi = 1d0 / dble(nevt(1)) / BMT
//C$$$               if(nevt(30).ne.0) SCALEf = 1d0 / dble(nevt(30)) / BMT
//C$$$         do 1016 iy=2,3
//C$$$           yra = 1d0
//C$$$           if(iy .eq. 2)yra = 2d0
//C$$$         do 1015 i=1,50
//C$$$           WRITE(nmfile(iy),251) BMT*dble(I - 0.5),
//C$$$     &     SCALEi*DMYil(I,iy)/yra, SCALEf*DMYfl(I,iy)/yra,
//C$$$     &     SCALEi*DMYis(I,iy)/yra, SCALEf*DMYfs(I,iy)/yra
//C$$$ 1015   continue
//C$$$ 1016 continue
//C$$$csp07/05 end
//C$$$clin-3/30/00 event-by-event average & variance of v2:
//C$$$           if(nevt(30).ge.1) then
//C$$$              do 1017 iy=1,3
//C$$$                 v2pavg(iy)=v2psum(iy)/nevt(30)
//C$$$                 v2var0=v2p2sm(iy)/nevt(30)-v2pavg(iy)**2
//C$$$clin-5/2012:
//C$$$c                 if(v2var0.gt.0) varv2p(iy)=dsqrt(v2var0)
//C$$$                 if(v2var0.gt.0d0) varv2p(iy)=dsqrt(v2var0)
//C$$$ 1017 continue
//C$$$              write(49, 240) 'EBE v2p,v2p(y2),v2p(y1): avg=', v2pavg
//C$$$              write(49, 240) 'EBE v2p,v2p(y2),v2p(y1): var=', varv2p
//C$$$           endif
//C$$$clin-8/2015:
//C$$$clin-11/28/00 for animation:
//C$$$           if(ifanim.eq.1) then
//C$$$              do 1018 I=1,mul
//C$$$                 if(FT5(I).le.t2time) then
//C$$$                    write(15,140) ITYP5(i),GX5(i),GY5(i),GZ5(i),FT5(i)
//C$$$                 endif
//C$$$ 1018         continue
//C$$$ 140          format(i10,4(2x,f7.2))
//C$$$clin-11/29/00 signal the end of animation file:
//C$$$              write(10,*) -10.
//C$$$              write(10,*) 0
//C$$$              write(11,*) -10.
//C$$$              write(11,*) 0
//C$$$              close(10)
//C$$$              close(11)
//C$$$              close(15)
//C$$$           endif
//C$$$
//C$$$clin-5/18/01 calculate v2 for active partons:
//C$$$           do 450 ianp=1,30
//C$$$              do 400 iy=1,3
//C$$$                 v2pact=0d0
//C$$$                 v2p2ac=0d0
//C$$$                 s2pact=0d0
//C$$$                 etpact=0d0
//C$$$                 xnacti=0d0
//C$$$clin-5/2012:
//C$$$c                 if(xnpf(ianp,iy).gt.1) then
//C$$$                 if(xnpf(ianp,iy).gt.1d0) then
//C$$$c     reconstruct the sum of v2p, v2p2, s2p, etp, and xnp for formed partons:
//C$$$                    v2pact=v2pf(ianp,iy)*xnpf(ianp,iy)
//C$$$                    v2p2ac=(v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)
//C$$$     1 +v2pf(ianp,iy)**2)*xnpf(ianp,iy)
//C$$$                    s2pact=s2pf(ianp,iy)*xnpf(ianp,iy)
//C$$$                    etpact=etpf(ianp,iy)*dble(nevt(ianp))
//C$$$                    xnpact=xnpf(ianp,iy)
//C$$$c
//C$$$                    do 1019 kanp=1,ianp
//C$$$                       v2pact=v2pact-v2pfrz(kanp,iy)
//C$$$                       v2p2ac=v2p2ac-v2p2fz(kanp,iy)
//C$$$                       s2pact=s2pact-s2pfrz(kanp,iy)
//C$$$                       etpact=etpact-etpfrz(kanp,iy)
//C$$$                       xnpact=xnpact-xnpfrz(kanp,iy)
//C$$$ 1019               continue
//C$$$c     save the sum of v2p, v2p2, s2p, etp, and xnp for formed partons:
//C$$$                    v2ph=v2pact
//C$$$                    v2ph2=v2p2ac
//C$$$                    s2ph=s2pact
//C$$$                    etph=etpact
//C$$$                    xnp2=xnpact/2d0
//C$$$c
//C$$$clin-5/2012:
//C$$$c                    if(xnpact.gt.1.and.nevt(ianp).ne.0) then
//C$$$                    if(xnpact.gt.1d0.and.nevt(ianp).ne.0) then
//C$$$                       v2pact=v2pact/xnpact
//C$$$                       v2p2ac=dsqrt((v2p2ac/xnpact
//C$$$     1                    -v2pact**2)/(xnpact-1))
//C$$$                       s2pact=s2pact/xnpact
//C$$$                       xnacti=dble(xnpact/dble(nevt(ianp)))
//C$$$                       etpact=etpact/dble(nevt(ianp))
//C$$$                       write (nfile(iy)+2, 250) tsp(ianp),v2pact,
//C$$$     1 v2p2ac,s2pact,etpact,xnacti,
//C$$$     2 xnpf(ianp,iy)/dble(nevt(ianp)),nevt(ianp)
//C$$$                    endif
//C$$$                 endif
//C$$$c     To calculate combined v2 for active partons plus formed hadrons,
//C$$$c     add the sum of v2h, v2h2, s2h, eth, and xnh for formed hadrons:
//C$$$c     scale the hadron part in case nevt(ianp) != nevent:
//C$$$                 shadr=dble(nevt(ianp))/dble(nevent)
//C$$$                 ianh=ianp
//C$$$                 v2ph=v2ph+v2h(ianh,iy)*xnhadr(ianh,iy)*shadr
//C$$$                 v2ph2=v2ph2+(v2h2(ianh,iy)**2*(xnhadr(ianh,iy)-1)
//C$$$     1 +v2h(ianh,iy)**2)*xnhadr(ianh,iy)*shadr
//C$$$                 s2ph=s2ph+s2h(ianh,iy)*xnhadr(ianh,iy)*shadr
//C$$$                 etph=etph+eth(ianh,iy)*dble(nevent)*shadr
//C$$$                 xnph=xnpact+xnhadr(ianh,iy)*shadr
//C$$$                 xnp2h=xnp2+xnhadr(ianh,iy)*shadr
//C$$$clin-8/2015 to avoid IEEE_DIVIDE_BY_ZERO:
//C$$$cclin-5/2012:
//C$$$cc                 if(xnph.gt.1) then
//C$$$c                 if(xnph.gt.1d0) then
//C$$$                 if(xnph.gt.1d0.and.nevt(ianp).ne.0) then
//C$$$                    v2ph=v2ph/xnph
//C$$$                    v2ph2=dsqrt((v2ph2/xnph-v2ph**2)/(xnph-1))
//C$$$                    s2ph=s2ph/xnph
//C$$$                    etph=etph/dble(nevt(ianp))
//C$$$                    xnp2=xnp2/dble(nevt(ianp))
//C$$$                    xnp2h=xnp2h/dble(nevent)
//C$$$clin-8/2015
//C$$$c                    if(tsp(ianp).le.(ntmax*dt))
//C$$$                    if(tsp(ianp).le.dble(ntmax*dt))
//C$$$     1                    write (nfile(iy)+3, 250) tsp(ianp),v2ph,
//C$$$     2 v2ph2,s2ph,etph,xnp2h,xnp2,nevt(ianp)
//C$$$                 endif
//C$$$c
//C$$$ 400              continue
//C$$$ 450       continue
//C$$$           do 550 ianp=1,30
//C$$$              do 500 iy=1,3
//C$$$                 v2pact=0d0
//C$$$                 v2p2ac=0d0
//C$$$                 s2pact=0d0
//C$$$                 etpact=0d0
//C$$$                 xnacti=0d0
//C$$$c     reconstruct the sum of v2p, v2p2, s2p, etp, and xnp for formed partons:
//C$$$                    v2pact=v2pf(ianp,iy)*xnpf(ianp,iy)
//C$$$                    v2p2ac=(v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)
//C$$$     1 +v2pf(ianp,iy)**2)*xnpf(ianp,iy)
//C$$$                    s2pact=s2pf(ianp,iy)*xnpf(ianp,iy)
//C$$$                    etpact=etpf(ianp,iy)*dble(nevt(ianp))
//C$$$                    xnpact=xnpf(ianp,iy)
//C$$$ 500              continue
//C$$$ 550           continue
//C$$$           close (620)
//C$$$           close (630)
//C$$$           do 1021 nf=1,3
//C$$$              do 1020 ifile=0,3
//C$$$                 close(nfile(nf)+ifile)
//C$$$ 1020        continue
//C$$$ 1021     continue
//C$$$           do 1022 nf=1,3
//C$$$              close(740+nf)
//C$$$ 1022      continue
//C$$$        endif
//C$$$ 200        format(2x,f5.2,3(2x,f7.4),2(2x,f9.2),i6,2x,f9.2)
//C$$$ 210        format(2x,f5.2,3(2x,f7.4),3(2x,f9.2))
//C$$$ 240        format(a30,3(2x,f9.5))
//C$$$ 250        format(2x,f5.2,3(2x,f7.4),3(2x,f9.2),i6)
//C$$$csp
//C$$$ 201        format(2x,f5.2,4(2x,f9.2))
//C$$$ 251        format(5e15.5)
//C$$$c
//C$$$        return
//C$$$        end
//C$$$
//C=======================================================================
//C     Calculate flow from formed hadrons, called by art1e.f:
//C     Note: numbers in art not in double precision!
void
flowh(
  common& cmn,
  float const& ct)
{
  FEM_CMN_SVE(flowh);
  common_write write(cmn);
  arr_ref<double, 2> v2h(cmn.v2h, dimension(30, 3));
  arr_ref<double, 2> xnhadr(cmn.xnhadr, dimension(30, 3));
  arr_ref<double, 2> eth(cmn.eth, dimension(30, 3));
  arr_ref<double, 2> v2h2(cmn.v2h2, dimension(30, 3));
  arr_ref<double, 2> s2h(cmn.s2h, dimension(30, 3));
  arr_ref<double> v2hp(cmn.v2hp, dimension(3));
  arr_ref<double> xnhadp(cmn.xnhadp, dimension(3));
  arr_ref<double> v2hsum(cmn.v2hsum, dimension(3));
  arr_ref<double> v2h2sm(cmn.v2h2sm, dimension(3));
  int& itimeh = cmn.itimeh;
  int& num = cmn.num;
  const int maxstr = 150001;
  arr_cref<float, 2> r(cmn.r, dimension(3, maxstr));
  arr_cref<float, 2> p(static_cast<common_bb&>(cmn).p, dimension(3, maxstr));
  arr_cref<float> e(cmn.e, dimension(maxstr));
  arr_cref<int> lb(cmn.lb, dimension(maxstr));
  const int maxr = 1;
  arr_cref<int> massr(cmn.massr, dim1(0, maxr));
  //
  float& ene = sve.ene;
  int& i = sve.i;
  int& ia = sve.ia;
  int& ianh = sve.ianh;
  int& ic = sve.ic;
  int& ifanim = sve.ifanim;
  int& ii = sve.ii;
  int& iloop = sve.iloop;
  int& iy = sve.iy;
  int& j = sve.j;
  int& mult = sve.mult;
  int& nhadrn = sve.nhadrn;
  double& pt2 = sve.pt2;
  float& px = sve.px;
  float& py = sve.py;
  float& rap = sve.rap;
  arr_ref<float> tsh(sve.tsh, dimension(31));
  double& v2hadr = sve.v2hadr;
  arr_ref<double> v2hevt(sve.v2hevt, dimension(3));
  float& xperp2 = sve.xperp2;
  //Cc      SAVE /hflow/
  //Cc      SAVE /ebe/
  //Cc      SAVE /lastt/
  //Cc      SAVE /RUN/
  //Cc      SAVE /AA/
  //Cc      SAVE /BB/
  //Cc      SAVE /CC/
  //Cc      SAVE /EE/
  //Cc      SAVE /RR/
  //Cc      SAVE /anim/
  //Cc      SAVE /AREVT/
  //C
  FEM_DO_SAFE(ii, 1, 31) {
    tsh(ii) = fem::ffloat(ii - 1);
  }
  //C
  FEM_DO_SAFE(ianh, 1, 30) {
    if ((ct + 0.0001f) < tsh(ianh + 1) && (ct + 0.0001f) >= tsh(ianh)) {
      if (ianh == itimeh) {
        goto statement_101;
      }
      ia = 0;
      FEM_DO_SAFE(j, 1, num) {
        mult = massr(j);
        ia += massr(j - 1);
        FEM_DO_SAFE(ic, 1, mult) {
          i = ia + ic;
          //C     5/04/01 exclude leptons and photons:
          if (fem::iabs(lb(i) - 10000) < 100) {
            goto statement_100;
          }
          px = p(1, i);
          py = p(2, i);
          pt2 = fem::pow2(fem::dble(px)) + fem::pow2(fem::dble(py));
          //C     2/18/00 Note: e(i) gives the mass in ART:
          ene = fem::sqrt(fem::pow2(e(i)) + fem::sngl(pt2) + fem::pow2(p(3,
            i)));
          rap = 0.5f * fem::alog((ene + p(3, i)) / (ene - p(3, i)));
          //Ctest off density with 2fm radius and z:(-0.1*t,0.1*t):
          //C                rt_now=sqrt(r(1,i)**2+r(2,i)**2)
          //C                gz_now=r(3,i)
          //C                zmax=0.1*ct
          //C                volume=3.1416*(2.**2)*2*zmax
          //C                if(rt_now.gt.2.or.abs(gz_now).gt.zmax)
          //C     1               goto 100
          iloop = 1;
          if (fem::abs(rap) <= 1) {
            iloop = 2;
            if (fem::abs(rap) <= 0.5f) {
              iloop = 3;
            }
          }
          FEM_DO_SAFE(iy, 1, iloop) {
            if (pt2 > 0e0) {
              v2hadr = (fem::pow2(fem::dble(px)) - fem::pow2(
                fem::dble(py))) / pt2;
              v2h(ianh, iy) += v2hadr;
              v2h2(ianh, iy) += fem::pow2(v2hadr);
              if (fem::dabs(v2hadr) > 1e0) {
                write(1, star), "v2hadr>1", v2hadr, px, py;
              }
            }
            xperp2 = fem::pow2(r(1, i)) + fem::pow2(r(2, i));
            if (xperp2 > 0.f) {
              s2h(ianh, iy) += fem::dble((fem::pow2(r(1, i)) - fem::pow2(r(2,
                i))) / xperp2);
            }
            eth(ianh, iy) += fem::dble(fem::sqrt(fem::pow2(e(i)) +
              fem::sngl(pt2)));
            //Ctest off density:
            //C               eth(ianh,iy)=eth(ianh,iy)
            //C     1                  +dble(SQRT(e(i)**2+sngl(pt2)+p(3,i)**2))/volume
            xnhadr(ianh, iy) += 1e0;
          }
          statement_100:;
        }
      }
      itimeh = ianh;
      //Clin-5/04/01 ebe v2 variables:
      if (ianh == 30) {
        FEM_DO_SAFE(iy, 1, 3) {
          nhadrn = idint(xnhadr(ianh, iy) - xnhadp(iy));
          if (nhadrn != 0) {
            v2hevt(iy) = (v2h(ianh, iy) - v2hp(iy)) / fem::dble(nhadrn);
            v2hsum(iy) += v2hevt(iy);
            v2h2sm(iy) += fem::pow2(v2hevt(iy));
            v2hp(iy) = v2h(ianh, iy);
            xnhadp(iy) = xnhadr(ianh, iy);
          }
        }
        write(88, "(i10,3(2x,f9.5))"), cmn.iaevt, v2hevt;
      }
      goto statement_101;
    }
  }
  //Clin-8/2015:
  //Clin-11/27/00 for animation:
  //Ctest off     turn off animation output (0 to turn off and 1 to turn on):
  statement_101:
  ifanim = 0;
  if (ifanim == 1) {
    ia = 0;
    FEM_DO_SAFE(j, 1, num) {
      mult = massr(j);
      ia += massr(j - 1);
      write(10, star), ct;
      write(10, star), mult;
      FEM_DO_SAFE(ic, 1, mult) {
        i = ia + ic;
        //Clin-6/2013 for animation:
        if (fem::amax1(fem::abs(r(1, i)), fem::abs(r(2, i)), fem::abs(r(3,
            i))) < 9999) {
          write(10, "(i6,7(1x,f9.3))"), lb(i), r(1, i), r(2, i), r(3,
            i), p(1, i), p(2, i), p(3, i), e(i);
        }
        else {
          write(10, "(i6,3(1x,e9.3),4(1x,f9.3))"), lb(i), r(1, i), r(2,
            i), r(3, i), p(1, i), p(2, i), p(3, i), e(i);
        }
      }
    }
    return;
  }
}

struct flowh0_save
{
  int ii;
  int iy;
  int nunit;
  arr<float> tsh;
  arr<double> v2havg;
  arr<double> varv2h;
  float xmulth;

  flowh0_save() :
    ii(fem::int0),
    iy(fem::int0),
    nunit(fem::int0),
    tsh(dimension(31), fem::fill0),
    v2havg(dimension(3), fem::fill0),
    varv2h(dimension(3), fem::fill0),
    xmulth(fem::float0)
  {}
};

//C
//C=======================================================================
void
flowh0(
  common& cmn,
  int const& nevnt,
  int const& idd)
{
  FEM_CMN_SVE(flowh0);
  common_write write(cmn);
  // COMMON hflow
  arr_ref<double, 2> v2h(cmn.v2h, dimension(30, 3));
  arr_ref<double, 2> xnhadr(cmn.xnhadr, dimension(30, 3));
  arr_ref<double, 2> eth(cmn.eth, dimension(30, 3));
  arr_ref<double, 2> v2h2(cmn.v2h2, dimension(30, 3));
  arr_ref<double, 2> s2h(cmn.s2h, dimension(30, 3));
  // COMMON ebe
  arr_ref<double> v2hp(cmn.v2hp, dimension(3));
  arr_ref<double> xnhadp(cmn.xnhadp, dimension(3));
  arr_ref<double> v2hsum(cmn.v2hsum, dimension(3));
  arr_ref<double> v2h2sm(cmn.v2h2sm, dimension(3));
  //
  // SAVE
  int& ii = sve.ii;
  int& iy = sve.iy;
  int& nunit = sve.nunit;
  arr_ref<float> tsh(sve.tsh, dimension(31));
  arr_ref<double> v2havg(sve.v2havg, dimension(3));
  arr_ref<double> varv2h(sve.varv2h, dimension(3));
  float& xmulth = sve.xmulth;
  //
  static const char* format_240 = "(a30,3(2x,f9.5))";
  //C
  //Cc      SAVE /hflow/
  //Cc      SAVE /ebe/
  //Cc      SAVE /input1/
  //Cc      SAVE /INPUT2/
  //Cc      SAVE /lastt/
  //C
  //C     idd=0: initialization for flow analysis, called by artdri.f::
  if (idd == 0) {
    cmn.itimeh = 0;
    //C
    FEM_DO_SAFE(ii, 1, 31) {
      tsh(ii) = fem::ffloat(ii - 1);
    }
    //C
    FEM_DO_SAFE(ii, 1, 30) {
      FEM_DO_SAFE(iy, 1, 3) {
        v2h(ii, iy) = 0e0;
        xnhadr(ii, iy) = 0e0;
        eth(ii, iy) = 0e0;
        v2h2(ii, iy) = 0e0;
        s2h(ii, iy) = 0e0;
      }
    }
    FEM_DO_SAFE(iy, 1, 3) {
      v2hp(iy) = 0e0;
      xnhadp(iy) = 0e0;
      v2hsum(iy) = 0e0;
      v2h2sm(iy) = 0e0;
      if (iy == 1) {
        nunit = 59;
      }
      else if (iy == 2) {
        nunit = 68;
      }
      else {
        nunit = 69;
      }
      write(nunit, star), "   tsh,   v2h,     v2h2,     s2h, " +
        str_cref(" eth,   xmulth");
    }
    //C     idd=2: calculate average hadronic elliptic flow, called by artdri.f:
  }
  else if (idd == 2) {
    FEM_DO_SAFE(ii, 1, 30) {
      FEM_DO_SAFE(iy, 1, 3) {
        if (xnhadr(ii, iy) == 0) {
          xmulth = 0.f;
        }
        else if (xnhadr(ii, iy) > 1) {
          v2h(ii, iy) = v2h(ii, iy) / xnhadr(ii, iy);
          eth(ii, iy) = eth(ii, iy) / fem::dble(nevnt);
          v2h2(ii, iy) = fem::dsqrt((v2h2(ii, iy) / xnhadr(ii, iy) -
            fem::pow2(v2h(ii, iy))) / (xnhadr(ii, iy) - 1));
          s2h(ii, iy) = s2h(ii, iy) / xnhadr(ii, iy);
          xmulth = fem::sngl(xnhadr(ii, iy) / nevnt);
        }
        if (iy == 1) {
          nunit = 59;
        }
        else if (iy == 2) {
          nunit = 68;
        }
        else {
          nunit = 69;
        }
        if (tsh(ii) <= (cmn.ntmax * cmn.dt)) {
          write(nunit, "(2x,f5.2,3(2x,f7.4),2(2x,f9.2))"), tsh(ii),
            v2h(ii, iy), v2h2(ii, iy), s2h(ii, iy), eth(ii, iy),
            xmulth;
        }
      }
    }
    //C     event-by-event average & variance of v2h:
    FEM_DO_SAFE(iy, 1, 3) {
      v2havg(iy) = v2hsum(iy) / fem::dble(nevnt);
      varv2h(iy) = fem::dsqrt(v2h2sm(iy) / fem::dble(nevnt) -
        fem::pow2(v2havg(iy)));
    }
    write(88, format_240), "EBE v2h,v2h(y2),v2h(y1): avg=", v2havg;
    write(88, format_240), "EBE v2h,v2h(y2),v2h(y1): var=", varv2h;
  }
}

struct iniflw_save
{
  int i;
  int ityp;
  int j;
  float pt2;
  float px;
  float py;
  float xh;
  float xm;
  float xt2;
  float yh;

  iniflw_save() :
    i(fem::int0),
    ityp(fem::int0),
    j(fem::int0),
    pt2(fem::float0),
    px(fem::float0),
    py(fem::float0),
    xh(fem::float0),
    xm(fem::float0),
    xt2(fem::float0),
    yh(fem::float0)
  {}
};

//C
//C=======================================================================
//C     2/23/00 flow from all initial hadrons just before entering ARTMN:
void
iniflw(
  common& cmn,
  int const& nevnt,
  int const& idd)
{
  FEM_CMN_SVE(iniflw);
  const int maxr = 1;
  arr_cref<int> multi1(cmn.multi1, dimension(maxr));
  const int maxstr = 150001;
  arr_cref<int, 2> ityp1(cmn.ityp1, dimension(maxstr, maxr));
  arr_cref<float, 2> gx1(cmn.gx1, dimension(maxstr, maxr));
  arr_cref<float, 2> gy1(cmn.gy1, dimension(maxstr, maxr));
  arr_cref<float, 2> px1(static_cast<common_arprc1&>(cmn).px1,
    dimension(maxstr, maxr));
  arr_cref<float, 2> py1(static_cast<common_arprc1&>(cmn).py1,
    dimension(maxstr, maxr));
  arr_cref<float, 2> xm1(cmn.xm1, dimension(maxstr, maxr));
  double& v2i = cmn.v2i;
  double& eti = cmn.eti;
  double& xmulti = cmn.xmulti;
  double& v2mi = cmn.v2mi;
  double& s2mi = cmn.s2mi;
  double& xmmult = cmn.xmmult;
  double& v2bi = cmn.v2bi;
  double& s2bi = cmn.s2bi;
  double& xbmult = cmn.xbmult;
  //
  int& i = sve.i;
  int& ityp = sve.ityp;
  int& j = sve.j;
  float& pt2 = sve.pt2;
  float& px = sve.px;
  float& py = sve.py;
  float& xh = sve.xh;
  float& xm = sve.xm;
  float& xt2 = sve.xt2;
  float& yh = sve.yh;
  //Cc      SAVE /RUN/
  //Cc      SAVE /ARERC1/
  //Cc      SAVE /ARPRC1/
  //Cc      SAVE /iflow/
  //C
  if (idd == 0) {
    v2i = 0e0;
    eti = 0e0;
    xmulti = 0e0;
    v2mi = 0e0;
    s2mi = 0e0;
    xmmult = 0e0;
    v2bi = 0e0;
    s2bi = 0e0;
    xbmult = 0e0;
  }
  else if (idd == 1) {
    FEM_DO_SAFE(j, 1, cmn.num) {
      FEM_DO_SAFE(i, 1, multi1(j)) {
        ityp = ityp1(i, j);
        //C     all hadrons:
        if (ityp >  - 100 && ityp < 100) {
          goto statement_100;
        }
        xmulti += 1.e0;
        px = px1(i, j);
        py = py1(i, j);
        xm = xm1(i, j);
        pt2 = fem::pow2(px) + fem::pow2(py);
        xh = gx1(i, j);
        yh = gy1(i, j);
        xt2 = fem::pow2(xh) + fem::pow2(yh);
        if (pt2 > 0) {
          v2i += fem::dble((fem::pow2(px) - fem::pow2(py)) / pt2);
        }
        eti += fem::dble(fem::sqrt(fem::pow2(px) + fem::pow2(py) +
          fem::pow2(xm)));
        //C     baryons only:
        if (ityp <  - 1000 || ityp > 1000) {
          xbmult += 1.e0;
          if (pt2 > 0) {
            v2bi += fem::dble((fem::pow2(px) - fem::pow2(py)) / pt2);
          }
          if (xt2 > 0) {
            s2bi += fem::dble((fem::pow2(xh) - fem::pow2(yh)) / xt2);
          }
          //C     mesons only:
        }
        else {
          xmmult += 1.e0;
          if (pt2 > 0) {
            v2mi += fem::dble((fem::pow2(px) - fem::pow2(py)) / pt2);
          }
          if (xt2 > 0) {
            s2mi += fem::dble((fem::pow2(xh) - fem::pow2(yh)) / xt2);
          }
        }
        statement_100:;
      }
    }
  }
  else if (idd == 2) {
    if (xmulti != 0) {
      v2i = v2i / xmulti;
    }
    eti = eti / fem::dble(nevnt);
    xmulti = xmulti / fem::dble(nevnt);
    if (xmmult != 0) {
      v2mi = v2mi / xmmult;
      s2mi = s2mi / xmmult;
    }
    xmmult = xmmult / fem::dble(nevnt);
    if (xbmult != 0) {
      v2bi = v2bi / xbmult;
      s2bi = s2bi / xbmult;
    }
    xbmult = xbmult / fem::dble(nevnt);
  }
  //C
}

struct frztm_save
{
  double detf;
  double detp;
  double dxnf;
  double dxnp;
  double etf;
  double eth0;
  double eth2;
  double etp;
  int ii;
  int ip;
  arr<double> tsf;
  double xnf;
  double xnp;

  frztm_save() :
    detf(fem::double0),
    detp(fem::double0),
    dxnf(fem::double0),
    dxnp(fem::double0),
    etf(fem::double0),
    eth0(fem::double0),
    eth2(fem::double0),
    etp(fem::double0),
    ii(fem::int0),
    ip(fem::int0),
    tsf(dimension(31), fem::fill0),
    xnf(fem::double0),
    xnp(fem::double0)
  {}
};

//C
//C=======================================================================
//C     2/25/00 dN/dt analysis for production (before ZPCMN)
//C     and freezeout (right after ZPCMN) for all partons.
void
frztm(
  common& cmn,
  int const& nevnt,
  int const& idd)
{
  FEM_CMN_SVE(frztm);
  common_write write(cmn);
  // COMMON prec1
  const int maxptn = 400001;
  arr_cref<double> ft0(cmn.ft0, dimension(maxptn));
  arr_cref<double> px0(cmn.px0, dimension(maxptn));
  arr_cref<double> py0(cmn.py0, dimension(maxptn));
  arr_cref<double> xmass0(cmn.xmass0, dimension(maxptn));
  // COMMON prec2
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_cref<double> px5(cmn.px5, dimension(maxptn));
  arr_cref<double> py5(cmn.py5, dimension(maxptn));
  arr_cref<double> xmass5(cmn.xmass5, dimension(maxptn));
  // COMMON frzout
  arr_ref<double> xnprod(cmn.xnprod, dimension(30));
  arr_ref<double> etprod(cmn.etprod, dimension(30));
  arr_ref<double> xnfrz(cmn.xnfrz, dimension(30));
  arr_ref<double> etfrz(cmn.etfrz, dimension(30));
  arr_ref<double> dnprod(cmn.dnprod, dimension(30));
  arr_ref<double> detpro(cmn.detpro, dimension(30));
  arr_ref<double> dnfrz(cmn.dnfrz, dimension(30));
  arr_ref<double> detfrz(cmn.detfrz, dimension(30));
  //
  // SAVE
  double& detf = sve.detf;
  double& detp = sve.detp;
  double& dxnf = sve.dxnf;
  double& dxnp = sve.dxnp;
  double& etf = sve.etf;
  double& eth0 = sve.eth0;
  double& eth2 = sve.eth2;
  double& etp = sve.etp;
  int& ii = sve.ii;
  int& ip = sve.ip;
  arr_ref<double> tsf(sve.tsf, dimension(31));
  double& xnf = sve.xnf;
  double& xnp = sve.xnp;
  //
  if (is_called_first_time) {
    {
      fem::data_values data;
      data.values, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f;
      data.values, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f;
      data.values, 1.6f, 1.7f, 1.8f, 1.9f, 2, 3, 4, 5;
      data.values, 6, 7, 8, 9, 10, 20, 30;
      data, tsf;
    }
  }
  static const char* format_200 = "(2x,f9.2,4(2x,f10.2))";
  //C
  //Cc      SAVE /PARA1/
  //Cc      SAVE /prec1/
  //Cc      SAVE /prec2/
  //Cc      SAVE /frzout/
  //C
  if (idd == 0) {
    FEM_DO_SAFE(ii, 1, 30) {
      xnprod(ii) = 0e0;
      etprod(ii) = 0e0;
      xnfrz(ii) = 0e0;
      etfrz(ii) = 0e0;
      dnprod(ii) = 0e0;
      detpro(ii) = 0e0;
      dnfrz(ii) = 0e0;
      detfrz(ii) = 0e0;
    }
    cmn.io.open(86, "ana1/production.dat")
      .status("UNKNOWN");
    cmn.io.open(87, "ana1/freezeout.dat")
      .status("UNKNOWN");
  }
  else if (idd == 1) {
    FEM_DO_SAFE(ip, 1, cmn.mul) {
      FEM_DO_SAFE(ii, 1, 30) {
        eth0 = fem::dsqrt(fem::pow2(px0(ip)) + fem::pow2(py0(ip)) +
          fem::pow2(xmass0(ip)));
        eth2 = fem::dsqrt(fem::pow2(px5(ip)) + fem::pow2(py5(ip)) +
          fem::pow2(xmass5(ip)));
        //C     total number and Et produced by time tsf(ii):
        if (ft0(ip) < tsf(ii + 1)) {
          xnprod(ii) += 1e0;
          etprod(ii) += eth0;
          //C     number and Et produced from time tsf(ii) to tsf(ii+1):
          if (ft0(ip) >= tsf(ii)) {
            dnprod(ii) += 1e0;
            detpro(ii) += eth0;
          }
        }
        //C     total number and Et freezed out by time tsf(ii):
        if (ft5(ip) < tsf(ii + 1)) {
          xnfrz(ii) += 1e0;
          etfrz(ii) += eth2;
          //C     number and Et freezed out from time tsf(ii) to tsf(ii+1):
          if (ft5(ip) >= tsf(ii)) {
            dnfrz(ii) += 1e0;
            detfrz(ii) += eth2;
          }
        }
      }
    }
  }
  else if (idd == 2) {
    write(86, star), "       t,       np,       dnp/dt,      etp " +
      str_cref(" detp/dt");
    write(87, star), "       t,       nf,       dnf/dt,      etf " +
      str_cref(" detf/dt");
    FEM_DO_SAFE(ii, 1, 30) {
      xnp = xnprod(ii) / fem::dble(nevnt);
      xnf = xnfrz(ii) / fem::dble(nevnt);
      etp = etprod(ii) / fem::dble(nevnt);
      etf = etfrz(ii) / fem::dble(nevnt);
      dxnp = dnprod(ii) / fem::dble(nevnt) / (tsf(ii + 1) - tsf(ii));
      dxnf = dnfrz(ii) / fem::dble(nevnt) / (tsf(ii + 1) - tsf(ii));
      detp = detpro(ii) / fem::dble(nevnt) / (tsf(ii + 1) - tsf(ii));
      detf = detfrz(ii) / fem::dble(nevnt) / (tsf(ii + 1) - tsf(ii));
      write(86, format_200), tsf(ii + 1), xnp, dxnp, etp, detp;
      write(87, format_200), tsf(ii + 1), xnf, dxnf, etf, detf;
    }
  }
  //C
}

struct minijet_out_save
{
  float ft;
  float gx;
  float gy;
  float gz;
  int i;
  int ityp;
  int j;
  float pt;
  float px;
  float py;
  float pz;
  float xmass;

  minijet_out_save() :
    ft(fem::float0),
    gx(fem::float0),
    gy(fem::float0),
    gz(fem::float0),
    i(fem::int0),
    ityp(fem::int0),
    j(fem::int0),
    pt(fem::float0),
    px(fem::float0),
    py(fem::float0),
    pz(fem::float0),
    xmass(fem::float0)
  {}
};

//C
//C=======================================================================
//Clin-6/2009 write out initial minijet information
//C     before propagating to its formation time:
//Clin-2/2012:
//C        subroutine minijet_out(BB)
void
minijet_out(
  common& cmn,
  float const& bb,
  float const& phirp)
{
  FEM_CMN_SVE(minijet_out);
  common_write write(cmn);
  // COMMON hparnt
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  // COMMON hjcrdn
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  // COMMON hjjet1
  arr_cref<int> npj(cmn.npj, dimension(300));
  arr_cref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_cref<float, 2> pjpx(cmn.pjpx, dimension(300, 500));
  arr_cref<float, 2> pjpy(cmn.pjpy, dimension(300, 500));
  arr_cref<float, 2> pjpz(cmn.pjpz, dimension(300, 500));
  arr_cref<float, 2> pjpm(cmn.pjpm, dimension(300, 500));
  arr_cref<int> ntj(cmn.ntj, dimension(300));
  arr_cref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  arr_cref<float, 2> pjtx(cmn.pjtx, dimension(300, 500));
  arr_cref<float, 2> pjty(cmn.pjty, dimension(300, 500));
  arr_cref<float, 2> pjtz(cmn.pjtz, dimension(300, 500));
  arr_cref<float, 2> pjtm(cmn.pjtm, dimension(300, 500));
  // COMMON hjjet2
  int& nsg = cmn.nsg;
  const int maxstr = 150001;
  arr_cref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_cref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_cref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_cref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_cref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_cref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_cref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  // COMMON para7
  int& ioscar = cmn.ioscar;
  // COMMON phidcy
  float& pttrig = cmn.pttrig;
  int& ntrig = cmn.ntrig;
  //
  // SAVE
  float& ft = sve.ft;
  float& gx = sve.gx;
  float& gy = sve.gy;
  float& gz = sve.gz;
  int& i = sve.i;
  int& ityp = sve.ityp;
  int& j = sve.j;
  float& pt = sve.pt;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& xmass = sve.xmass;
  //
  static const char* format_200 =
    "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,2(1x,f8.2),2(2x,f2.0),2x,i2)";
  static const char* format_201 =
    "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,2(1x,e8.2),2(2x,f2.0),2x,i2)";
  ntrig = 0;
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      pt = fem::sqrt(fem::pow2(pjpx(i, j)) + fem::pow2(pjpy(i, j)));
      if (pt >= pttrig) {
        ntrig++;
      }
    }
  }
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      pt = fem::sqrt(fem::pow2(pjtx(i, j)) + fem::pow2(pjty(i, j)));
      if (pt >= pttrig) {
        ntrig++;
      }
    }
  }
  FEM_DO_SAFE(i, 1, nsg) {
    FEM_DO_SAFE(j, 1, njsg(i)) {
      pt = fem::sqrt(fem::pow2(pxsg(i, j)) + fem::pow2(pysg(i, j)));
      if (pt >= pttrig) {
        ntrig++;
      }
    }
  }
  //C     Require at least 1 initial minijet parton above the trigger Pt value:
  if (ntrig == 0) {
    return;
  }
  //C
  //C.....transfer data from HIJING to ZPC
  if (ioscar == 3) {
    write(96, star), cmn.iaevt, cmn.miss, ihnt2(1), ihnt2(3);
  }
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      ityp = kfpj(i, j);
      //C     write out not only gluons:
      //C              if(ityp.ne.21) goto 1007
      //Clin-2/2012:
      //C              gx=YP(1,I)+0.5*BB
      //C              gy=YP(2,I)
      gx = yp(1, i) + 0.5f * bb * fem::cos(phirp);
      gy = yp(2, i) + 0.5f * bb * fem::sin(phirp);
      gz = 0.f;
      ft = 0.f;
      px = pjpx(i, j);
      py = pjpy(i, j);
      pz = pjpz(i, j);
      xmass = pjpm(i, j);
      if (ioscar == 3) {
        if (fem::amax1(fem::abs(gx), fem::abs(gy), fem::abs(gz),
            fem::abs(ft)) < 9999) {
          write(96, format_200), ityp, px, py, pz, xmass, gx, gy, gz, ft, 1;
        }
        else {
          write(96, format_201), ityp, px, py, pz, xmass, gx, gy, gz, ft, 1;
        }
      }
    }
  }
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      ityp = kftj(i, j);
      //C              if(ityp.ne.21) goto 1009
      //Clin-2/2012:
      //C              gx=YT(1,I)-0.5*BB
      //C              gy=YT(2,I)
      gx = yt(1, i) - 0.5f * bb * fem::cos(phirp);
      gy = yt(2, i) - 0.5f * bb * fem::sin(phirp);
      gz = 0.f;
      ft = 0.f;
      px = pjtx(i, j);
      py = pjty(i, j);
      pz = pjtz(i, j);
      xmass = pjtm(i, j);
      if (ioscar == 3) {
        if (fem::amax1(fem::abs(gx), fem::abs(gy), fem::abs(gz),
            fem::abs(ft)) < 9999) {
          write(96, format_200), ityp, px, py, pz, xmass, gx, gy, gz, ft, 2;
        }
        else {
          write(96, format_201), ityp, px, py, pz, xmass, gx, gy, gz, ft, 2;
        }
      }
    }
  }
  FEM_DO_SAFE(i, 1, nsg) {
    FEM_DO_SAFE(j, 1, njsg(i)) {
      ityp = k2sg(i, j);
      //C              if(ityp.ne.21) goto 1011
      gx = 0.5f * (yp(1, iasg(i, 1)) + yt(1, iasg(i, 2)));
      gy = 0.5f * (yp(2, iasg(i, 1)) + yt(2, iasg(i, 2)));
      gz = 0.f;
      ft = 0.f;
      px = pxsg(i, j);
      py = pysg(i, j);
      pz = pzsg(i, j);
      xmass = pmsg(i, j);
      if (ioscar == 3) {
        if (fem::amax1(fem::abs(gx), fem::abs(gy), fem::abs(gz),
            fem::abs(ft)) < 9999) {
          write(96, format_200), ityp, px, py, pz, xmass, gx, gy, gz, ft, 3;
        }
        else {
          write(96, format_201), ityp, px, py, pz, xmass, gx, gy, gz, ft, 3;
        }
      }
    }
  }
  //C
}

struct embedhighpt_save
{
  int idpi;
  int idpi1;
  int idpis;
  int idqembd;
  int idqsoft;
  int idsart;
  int ipion;
  int ispion;
  int ixy;
  float phi;
  float pimass;
  float ptpi;
  float ptq;
  float pxpi;
  float pxpi1;
  float pxspi;
  float pypi;
  float pypi1;
  float pyspi;
  float pzpi;
  float pzpi1;
  float pzspi;
  float theta;
  float xjet;
  float xmq;
  float xmqsoft;
  float yjet;

  embedhighpt_save() :
    idpi(fem::int0),
    idpi1(fem::int0),
    idpis(fem::int0),
    idqembd(fem::int0),
    idqsoft(fem::int0),
    idsart(fem::int0),
    ipion(fem::int0),
    ispion(fem::int0),
    ixy(fem::int0),
    phi(fem::float0),
    pimass(fem::float0),
    ptpi(fem::float0),
    ptq(fem::float0),
    pxpi(fem::float0),
    pxpi1(fem::float0),
    pxspi(fem::float0),
    pypi(fem::float0),
    pypi1(fem::float0),
    pyspi(fem::float0),
    pzpi(fem::float0),
    pzpi1(fem::float0),
    pzspi(fem::float0),
    theta(fem::float0),
    xjet(fem::float0),
    xmq(fem::float0),
    xmqsoft(fem::float0),
    yjet(fem::float0)
  {}
};

//C
//C=======================================================================
//Clin-6/2009 embed back-to-back high-Pt quark/antiquark pair
//C     via embedding back-to-back high-Pt pion pair then melting the pion pair
//C     by generating the internal quark and antiquark momentum parallel to
//C      the pion momentum (in order to produce a high-Pt and a low Pt parton):
void
embedhighpt(
  common& cmn)
{
  FEM_CMN_SVE(embedhighpt);
  common_read read(cmn);
  common_write write(cmn);
  // COMMON embed
  int& iembed = cmn.iembed;
  int& nsembd = cmn.nsembd;
  float& pxqembd = cmn.pxqembd;
  float& pyqembd = cmn.pyqembd;
  float& psembd = cmn.psembd;
  float& phidecomp = cmn.phidecomp;
  // COMMON rndf77
  int& nseed = cmn.nseed;
  // COMMON hmain1
  float& eatt = cmn.eatt;
  int& natt = cmn.natt;
  // COMMON hmain2
  const int maxstr = 150001;
  arr_ref<int, 2> katt(cmn.katt, dimension(maxstr, 4));
  arr_ref<float, 2> patt(cmn.patt, dimension(maxstr, 4));
  // COMMON arprc
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
  // COMMON xyembed
  int& nxyjet = cmn.nxyjet;
  const int nxymax = 10001;
  arr_cref<float, 2> xyjet(cmn.xyjet, dimension(nxymax, 2));
  //
  // SAVE
  int& idpi = sve.idpi;
  int& idpi1 = sve.idpi1;
  int& idpis = sve.idpis;
  int& idqembd = sve.idqembd;
  int& idqsoft = sve.idqsoft;
  int& idsart = sve.idsart;
  int& ipion = sve.ipion;
  int& ixy = sve.ixy;
  float& phi = sve.phi;
  float& pimass = sve.pimass;
  float& ptpi = sve.ptpi;
  float& ptq = sve.ptq;
  float& pxpi = sve.pxpi;
  float& pxpi1 = sve.pxpi1;
  float& pxspi = sve.pxspi;
  float& pypi = sve.pypi;
  float& pypi1 = sve.pypi1;
  float& pyspi = sve.pyspi;
  float& pzpi = sve.pzpi;
  float& pzpi1 = sve.pzpi1;
  float& pzspi = sve.pzspi;
  float& theta = sve.theta;
  float& xjet = sve.xjet;
  float& xmq = sve.xmq;
  float& xmqsoft = sve.xmqsoft;
  float& yjet = sve.yjet;
  //
  //C
  if (iembed == 1 || iembed == 2) {
    xjet = cmn.xembd;
    yjet = cmn.yembd;
  }
  else if (iembed == 3 || iembed == 4) {
    if (cmn.nevent <= nxyjet) {
      read(97, star), xjet, yjet;
    }
    else {
      ixy = fem::mod(cmn.iaevt, nxyjet);
      if (ixy == 0) {
        ixy = nxyjet;
      }
      xjet = xyjet(ixy, 1);
      yjet = xyjet(ixy, 2);
    }
  }
  else {
    return;
  }
  //C
  ptq = fem::sqrt(fem::pow2(pxqembd) + fem::pow2(pyqembd));
  const float pichmass = 0.140f;
  if (ptq < (pichmass / 2.f)) {
    write(6, star), "Embedded quark transverse momentum is too small";
    FEM_STOP(0);
  }
  //C     Randomly embed u/ubar or d/dbar at high Pt:
  idqembd = 1 + fem::fint(2 * ranart(nseed));
  //C     Flavor content for the charged pion that contains the leading quark:
  if (idqembd == 1) {
    idqsoft = -2;
    idpi1 = -211;
  }
  else if (idqembd == 2) {
    idqsoft = -1;
    idpi1 = 211;
  }
  else {
    write(6, star), "Wrong quark flavor embedded";
    FEM_STOP(0);
  }
  //C     Caculate transverse momentum of the parent charged pion:
  xmq = ulmass(idqembd);
  xmqsoft = ulmass(idqsoft);
  ptpi = ((fem::pow2(pichmass) + fem::pow2(xmq) - fem::pow2(xmqsoft)) *
    ptq - fem::sqrt((fem::pow2(xmq) + fem::pow2(ptq)) * (fem::pow4(
    pichmass) - 2.f * fem::pow2(pichmass) * (fem::pow2(xmq) + fem::pow2(
    xmqsoft)) + fem::pow2((fem::pow2(xmq) - fem::pow2(xmqsoft)))))) / (
    2.f * fem::pow2(xmq));
  const float pi = 3.1415926f;
  if (iembed == 1 || iembed == 3) {
    pxpi1 = ptpi * pxqembd / ptq;
    pypi1 = ptpi * pyqembd / ptq;
    phidecomp = fem::acos(pxqembd / ptq);
    if (pyqembd < 0) {
      phidecomp = 2.f * pi - phidecomp;
    }
  }
  else {
    phidecomp = 2.f * pi * ranart(nseed);
    pxpi1 = ptpi * fem::cos(phidecomp);
    pypi1 = ptpi * fem::sin(phidecomp);
  }
  //C     Embedded quark/antiquark are assumed to have pz=0:
  pzpi1 = 0.f;
  //C     Insert the two parent charged pions,
  //C     ipion=1 for the pion containing the leading quark,
  //C     ipion=2 for the pion containing the leading antiquark of the same flavor:
  FEM_DO_SAFE(ipion, 1, 2) {
    if (ipion == 1) {
      idpi = idpi1;
      pxpi = pxpi1;
      pypi = pypi1;
      pzpi = pzpi1;
    }
    else if (ipion == 2) {
      idpi = -idpi1;
      pxpi = -pxpi1;
      pypi = -pypi1;
      pzpi = -pzpi1;
    }
    natt++;
    katt(natt, 1) = idpi;
    katt(natt, 2) = 40;
    katt(natt, 3) = 0;
    patt(natt, 1) = pxpi;
    patt(natt, 2) = pypi;
    patt(natt, 3) = pzpi;
    patt(natt, 4) = fem::sqrt(fem::pow2(pxpi) + fem::pow2(pypi) +
      fem::pow2(pzpi) + fem::pow2(pichmass));
    eatt += patt(natt, 4);
    gxar(natt) = xjet;
    gyar(natt) = yjet;
    gzar(natt) = 0.f;
    ftar(natt) = 0.f;
    itypar(natt) = katt(natt, 1);
    pxar(natt) = patt(natt, 1);
    pyar(natt) = patt(natt, 2);
    pzar(natt) = patt(natt, 3);
    pear(natt) = patt(natt, 4);
    xmar(natt) = pichmass;
  }
  //C
  //Clin-8/2009
  //C     Randomly embed a number of soft pions around each high-Pt quark in pair:
  const float pi0mass = 0.135f;
  if (nsembd > 0) {
    FEM_DO_SAFE(ipion, 1, 2) {
      FEM_DO_SAFE(sve.ispion, 1, nsembd) {
        idsart = 3 + fem::fint(3 * ranart(nseed));
        if (idsart == 3) {
          pimass = pichmass;
          idpis = -211;
        }
        else if (idsart == 4) {
          pimass = pi0mass;
          idpis = 111;
        }
        else {
          pimass = pichmass;
          idpis = 211;
        }
        natt++;
        katt(natt, 1) = idpis;
        katt(natt, 2) = 40;
        katt(natt, 3) = 0;
        //C     theta: relative angle between soft pion & associated high-Pt q or qbar,
        //C     generate theta and phi uniformly:
        //C     Note: it is not generated uniformly in solid angle because that gives
        //C     a valley at theta=0, unlike the jet-like correlation (a peak at theta=0).
        theta = cmn.tmaxembd * ranart(nseed);
        phi = 2.f * pi * ranart(nseed);
        pxspi = psembd * fem::sin(theta) * fem::cos(phi);
        pyspi = psembd * fem::sin(theta) * fem::sin(phi);
        pzspi = psembd * fem::cos(theta);
        if (ipion == 1) {
          rotate(pxpi1, pypi1, pzpi1, pxspi, pyspi, pzspi);
        }
        else {
          rotate(-pxpi1, -pypi1, -pzpi1, pxspi, pyspi, pzspi);
        }
        //Ctest off
        //C               write(99,*) "2  ",pxspi,pyspi,pzspi
        patt(natt, 1) = pxspi;
        patt(natt, 2) = pyspi;
        patt(natt, 3) = pzspi;
        patt(natt, 4) = fem::sqrt(fem::pow2(psembd) + fem::pow2(pimass));
        eatt += patt(natt, 4);
        gxar(natt) = xjet;
        gyar(natt) = yjet;
        gzar(natt) = 0.f;
        ftar(natt) = 0.f;
        itypar(natt) = katt(natt, 1);
        pxar(natt) = patt(natt, 1);
        pyar(natt) = patt(natt, 2);
        pzar(natt) = patt(natt, 3);
        pear(natt) = patt(natt, 4);
        xmar(natt) = pimass;
      }
    }
  }
  //Clin-8/2009-end
  //C
}

} // namespace AMPT
