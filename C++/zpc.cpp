#include <fem.hpp> // Fortran EMulation library of fable module

namespace AMPT {

using namespace fem::major_types;

void
inifrz(...)
{
  throw std::runtime_error(
    "Missing function implementation: inifrz");
}

void
local(...)
{
  throw std::runtime_error(
    "Missing function implementation: local");
}

void
zpstrg(...)
{
  throw std::runtime_error(
    "Missing function implementation: zpstrg");
}

struct common_para1
{
  int mul;

  common_para1() :
    mul(fem::int0)
  {}
};

struct common_para3
{
  int nsevt;
  int nevnt;
  int nsbrun;
  int ievt;
  int isbrun;

  common_para3() :
    nsevt(fem::int0),
    nevnt(fem::int0),
    nsbrun(fem::int0),
    ievt(fem::int0),
    isbrun(fem::int0)
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

struct common_rndm1
{
  int number;

  common_rndm1() :
    number(fem::int0)
  {}
};

struct common_rndm3
{
  int iseedp;

  common_rndm3() :
    iseedp(fem::int0)
  {}
};

struct common_ilist3
{
  double size1;
  double size2;
  double size3;
  double v1;
  double v2;
  double v3;
  double size;

  common_ilist3() :
    size1(fem::double0),
    size2(fem::double0),
    size3(fem::double0),
    v1(fem::double0),
    v2(fem::double0),
    v3(fem::double0),
    size(fem::double0)
  {}
};

struct common_para2
{
  double xmp;
  double xmu;
  double alpha;
  double rscut2;
  double cutof2;

  common_para2() :
    xmp(fem::double0),
    xmu(fem::double0),
    alpha(fem::double0),
    rscut2(fem::double0),
    cutof2(fem::double0)
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

struct common_para5
{
  int iconfg;
  int iordsc;

  common_para5() :
    iconfg(fem::int0),
    iordsc(fem::int0)
  {}
};

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

struct common_para6
{
  double centy;

  common_para6() :
    centy(fem::double0)
  {}
};

struct common_para4
{
  int iftflg;
  int ireflg;
  int igeflg;
  int ibstfg;

  common_para4() :
    iftflg(fem::int0),
    ireflg(fem::int0),
    igeflg(fem::int0),
    ibstfg(fem::int0)
  {}
};

struct common_par1
{
  double formt;

  common_par1() :
    formt(fem::double0)
  {}
};

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

struct common_ilist4
{
  static const int maxptn = 400001;

  int ifmpt;
  int ichkpt;
  arr<int> indx;

  common_ilist4() :
    ifmpt(fem::int0),
    ichkpt(fem::int0),
    indx(dimension(maxptn), fem::fill0)
  {}
};

const int common_ilist4::maxptn;

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

struct common_prec2
{
  static const int maxptn = 400001;

  arr<double> gx;
  arr<double> gy;
  arr<double> gz;
  arr<double> ft;
  arr<double> px;
  arr<double> py;
  arr<double> pz;
  arr<double> e;
  arr<double> xmass;
  arr<int> ityp;

  common_prec2() :
    gx(dimension(maxptn), fem::fill0),
    gy(dimension(maxptn), fem::fill0),
    gz(dimension(maxptn), fem::fill0),
    ft(dimension(maxptn), fem::fill0),
    px(dimension(maxptn), fem::fill0),
    py(dimension(maxptn), fem::fill0),
    pz(dimension(maxptn), fem::fill0),
    e(dimension(maxptn), fem::fill0),
    xmass(dimension(maxptn), fem::fill0),
    ityp(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec2::maxptn;

struct common_prec3
{
  static const int maxptn = 400001;

  arr<double> gxs;
  arr<double> gys;
  arr<double> gzs;
  arr<double> fts;
  arr<double> pxs;
  arr<double> pys;
  arr<double> pzs;
  arr<double> es;
  arr<double> xmasss;
  arr<int> ityps;

  common_prec3() :
    gxs(dimension(maxptn), fem::fill0),
    gys(dimension(maxptn), fem::fill0),
    gzs(dimension(maxptn), fem::fill0),
    fts(dimension(maxptn), fem::fill0),
    pxs(dimension(maxptn), fem::fill0),
    pys(dimension(maxptn), fem::fill0),
    pzs(dimension(maxptn), fem::fill0),
    es(dimension(maxptn), fem::fill0),
    xmasss(dimension(maxptn), fem::fill0),
    ityps(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec3::maxptn;

struct common_prec6
{
  static const int maxptn = 400001;

  arr<double> etas;
  arr<double> raps;
  arr<double> taus;

  common_prec6() :
    etas(dimension(maxptn), fem::fill0),
    raps(dimension(maxptn), fem::fill0),
    taus(dimension(maxptn), fem::fill0)
  {}
};

const int common_prec6::maxptn;

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

struct common_smearz
{
  double smearp;
  double smearh;

  common_smearz() :
    smearp(fem::double0),
    smearh(fem::double0)
  {}
};

struct common_precpb
{
  static const int maxptn = 400001;

  arr<double> vxp;
  arr<double> vyp;
  arr<double> vzp;

  common_precpb() :
    vxp(dimension(maxptn), fem::fill0),
    vyp(dimension(maxptn), fem::fill0),
    vzp(dimension(maxptn), fem::fill0)
  {}
};

const int common_precpb::maxptn;

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

struct common_ilist1
{
  static const int maxptn = 400001;

  int iscat;
  int jscat;
  arr<int> next;
  arr<int> last;
  int ictype;
  arr<int> icsta;
  arr<int> nic;
  arr<int> icels;

  common_ilist1() :
    iscat(fem::int0),
    jscat(fem::int0),
    next(dimension(maxptn), fem::fill0),
    last(dimension(maxptn), fem::fill0),
    ictype(fem::int0),
    icsta(dimension(maxptn), fem::fill0),
    nic(dimension(maxptn), fem::fill0),
    icels(dimension(maxptn), fem::fill0)
  {}
};

const int common_ilist1::maxptn;

struct common_ilist2
{
  int icell;
  arr<int, 3> icel;

  common_ilist2() :
    icell(fem::int0),
    icel(dimension(10, 10, 10), fem::fill0)
  {}
};

struct common_ilist6
{
  double t;
  int iopern;
  int icolln;

  common_ilist6() :
    t(fem::double0),
    iopern(fem::int0),
    icolln(fem::int0)
  {}
};

struct common_ana2
{
  arr<double> det;
  arr<double> dn;
  arr<double> detdy;
  arr<double> detdn;
  arr<double> dndy;
  arr<double> det1;
  arr<double> dn1;
  arr<double> detdy1;
  arr<double> detdn1;
  arr<double> dndy1;
  arr<double> det2;
  arr<double> dn2;
  arr<double> detdy2;
  arr<double> detdn2;
  arr<double> dndy2;

  common_ana2() :
    det(dimension(12), fem::fill0),
    dn(dimension(12), fem::fill0),
    detdy(dimension(12), fem::fill0),
    detdn(dimension(12), fem::fill0),
    dndy(dimension(12), fem::fill0),
    det1(dimension(12), fem::fill0),
    dn1(dimension(12), fem::fill0),
    detdy1(dimension(12), fem::fill0),
    detdn1(dimension(12), fem::fill0),
    dndy1(dimension(12), fem::fill0),
    det2(dimension(12), fem::fill0),
    dn2(dimension(12), fem::fill0),
    detdy2(dimension(12), fem::fill0),
    detdn2(dimension(12), fem::fill0),
    dndy2(dimension(12), fem::fill0)
  {}
};

struct common_aurec1
{
  int jxa;
  int jya;
  int jza;

  common_aurec1() :
    jxa(fem::int0),
    jya(fem::int0),
    jza(fem::int0)
  {}
};

struct common_aurec2
{
  static const int maxptn = 400001;

  arr<double> dgxa;
  arr<double> dgya;
  arr<double> dgza;

  common_aurec2() :
    dgxa(dimension(maxptn), fem::fill0),
    dgya(dimension(maxptn), fem::fill0),
    dgza(dimension(maxptn), fem::fill0)
  {}
};

const int common_aurec2::maxptn;

struct common_rndm2
{
  int iff;

  common_rndm2() :
    iff(fem::int0)
  {}
};

struct common_ana1
{
  arr<double> ts;

  common_ana1() :
    ts(dimension(12), fem::fill0)
  {}
};

struct common_ana3
{
  arr<double, 3> em;

  common_ana3() :
    em(dimension(4, 4, 12), fem::fill0)
  {}
};

struct common_ana4
{
  arr<double> fdetdy;
  arr<double> fdndy;
  arr<double> fdndpt;

  common_ana4() :
    fdetdy(dimension(24), fem::fill0),
    fdndy(dimension(24), fem::fill0),
    fdndpt(dimension(12), fem::fill0)
  {}
};

struct common :
  fem::common,
  common_para1,
  common_para3,
  common_prec1,
  common_rndm1,
  common_rndm3,
  common_ilist3,
  common_para2,
  common_lor,
  common_para5,
  common_prec5,
  common_para6,
  common_para4,
  common_par1,
  common_prec4,
  common_ilist4,
  common_ilist5,
  common_anim,
  common_prec2,
  common_prec3,
  common_prec6,
  common_ilist7,
  common_ilist8,
  common_smearz,
  common_precpb,
  common_precpa,
  common_frzprc,
  common_para7,
  common_arevt,
  common_ilist1,
  common_ilist2,
  common_ilist6,
  common_ana2,
  common_aurec1,
  common_aurec2,
  common_rndm2,
  common_ana1,
  common_ana3,
  common_ana4
{
  fem::variant_core common_cprod;
  fem::cmn_sve readi_sve;
  fem::cmn_sve ran1_sve;
  fem::cmn_sve posit1_sve;
  fem::cmn_sve posit2_sve;
  fem::cmn_sve posit3_sve;
  fem::cmn_sve energy_sve;
  fem::cmn_sve momntm_sve;
  fem::cmn_sve lorenz_sve;
  fem::cmn_sve genei_sve;
  fem::cmn_sve boosti_sve;
  fem::cmn_sve index1_sve;
  fem::cmn_sve ftime1_sve;
  fem::cmn_sve ftime_sve;
  fem::cmn_sve inirec_sve;
  fem::cmn_sve iilist_sve;
  fem::cmn_sve inian2_sve;
  fem::cmn_sve getict_sve;
  fem::cmn_sve newcre_sve;
  fem::cmn_sve celasn_sve;
  fem::cmn_sve oldcre_sve;
  fem::cmn_sve wallc1_sve;
  fem::cmn_sve wallc2_sve;
  fem::cmn_sve wallcb_sve;
  fem::cmn_sve fixtim_sve;
  fem::cmn_sve isco1_sve;
  fem::cmn_sve isco2_sve;
  fem::cmn_sve isco3_sve;
  fem::cmn_sve isco4_sve;
  fem::cmn_sve isco5_sve;
  fem::cmn_sve isco6_sve;
  fem::cmn_sve isco7_sve;
  fem::cmn_sve isco8_sve;
  fem::cmn_sve isco9_sve;
  fem::cmn_sve isco10_sve;
  fem::cmn_sve isco11_sve;
  fem::cmn_sve isco12_sve;
  fem::cmn_sve isco_sve;
  fem::cmn_sve mintm_sve;
  fem::cmn_sve chcell_sve;
  fem::cmn_sve chout_sve;
  fem::cmn_sve chin1_sve;
  fem::cmn_sve chin2_sve;
  fem::cmn_sve chin3_sve;
  fem::cmn_sve reor_sve;
  fem::cmn_sve dchcel_sve;
  fem::cmn_sve dchout_sve;
  fem::cmn_sve dchin1_sve;
  fem::cmn_sve dchin2_sve;
  fem::cmn_sve dchin3_sve;
  fem::cmn_sve cellre_sve;
  fem::cmn_sve newpos_sve;
  fem::cmn_sve getht_sve;
  fem::cmn_sve cropro_sve;
  fem::cmn_sve xnormv_sve;
  fem::cmn_sve zprota_sve;
  fem::cmn_sve newmom_sve;
  fem::cmn_sve ud2_sve;
  fem::cmn_sve chkcel_sve;
  fem::cmn_sve chkout_sve;
  fem::cmn_sve chkin1_sve;
  fem::cmn_sve chkin2_sve;
  fem::cmn_sve chkin3_sve;
  fem::cmn_sve ulist1_sve;
  fem::cmn_sve ulist_sve;
  fem::cmn_sve zpcrun_sve;
  fem::cmn_sve zpca1c_sve;
  fem::cmn_sve zpca1a_sve;
  fem::cmn_sve zpca2a_sve;
  fem::cmn_sve zpca2b_sve;
  fem::cmn_sve zpca2c_sve;
  fem::cmn_sve zpcou1_sve;
  fem::cmn_sve zpcou2_sve;
  fem::cmn_sve zpcmn_sve;
  fem::cmn_sve blockdata_zpcbdt_sve;
  fem::cmn_sve readpa_sve;
  fem::cmn_sve inian1_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct readi_save
{
  arr<double> field;
  int i;
  int neve;
  int ntyp;

  readi_save() :
    field(dimension(9), fem::fill0),
    i(fem::int0),
    neve(fem::int0),
    ntyp(fem::int0)
  {}
};

void
readi(
  common& cmn)
{
  FEM_CMN_SVE(readi);
  common_read read(cmn);
  int& mul = cmn.mul;
  int& nsevt = cmn.nsevt;
  int& ievt = cmn.ievt;
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
  //
  arr_ref<double> field(sve.field, dimension(9));
  int& i = sve.i;
  int& neve = sve.neve;
  int& ntyp = sve.ntyp;
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para3/
  //Cc      SAVE /prec1/
  FEM_DO_SAFE(i, 1, maxptn) {
    if (ievt != 1 && i == 1) {
      ityp0(i) = ntyp;
      gx0(1) = field(1);
      gy0(1) = field(2);
      gz0(1) = field(3);
      ft0(1) = field(4);
      px0(1) = field(5);
      py0(1) = field(6);
      pz0(1) = field(7);
      e0(1) = field(8);
      xmass0(i) = field(9);
      mul = 1;
    }
    else {
      statement_900:
      try {
        read(27, star), neve, ntyp, field;
      }
      catch (fem::read_end const&) {
        goto statement_1000;
      }
      if (neve < nsevt) {
        goto statement_900;
      }
      if (neve > nsevt + ievt - 1) {
        goto statement_1000;
      }
      ityp0(i) = ntyp;
      gx0(i) = field(1);
      gy0(i) = field(2);
      gz0(i) = field(3);
      ft0(i) = field(4);
      px0(i) = field(5);
      py0(i) = field(6);
      pz0(i) = field(7);
      e0(i) = field(8);
      xmass0(i) = field(9);
      mul++;
    }
  }
  //C
  statement_1000:;
  //C
}

struct ran1_save
{
  int iff;
  int ix1;
  int ix2;
  int ix3;
  int j;
  arr<double> r;

  ran1_save() :
    iff(fem::int0),
    ix1(fem::int0),
    ix2(fem::int0),
    ix3(fem::int0),
    j(fem::int0),
    r(dimension(97), fem::fill0)
  {}
};

double
ran1(
  common& cmn,
  int& idum)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(ran1);
  common_write write(cmn);
  // COMMON rndm1
  int& number = cmn.number;
  //
  // SAVE
  int& iff = sve.iff;
  int& ix1 = sve.ix1;
  int& ix2 = sve.ix2;
  int& ix3 = sve.ix3;
  int& j = sve.j;
  arr_ref<double> r(sve.r, dimension(97));
  //
  if (is_called_first_time) {
    iff = 0;
  }
  //C
  //C     return a uniform random deviate between 0.0 and 1.0. set idum to
  //C     any negative value to initialize or reinitialize the sequence.
  //C
  //Cc      SAVE /rndm1/
  //Clin-6/23/00 save ix1-3:
  //Clin-10/30/02 r unsaved, causing wrong values for ran1 when compiled with f77:
  //Cc      SAVE ix1,ix2,ix3,r
  //C
  const int ic1 = 54773;
  const int m1 = 259200;
  const int ia1 = 7141;
  const int m2 = 134456;
  const int m3 = 243000;
  const int ia2 = 8121;
  const int ic2 = 28411;
  const double rm2 = 1e0 / m2;
  const double rm1 = 1e0 / m1;
  if (idum < 0 || iff == 0) {
    iff = 1;
    ix1 = fem::mod(ic1 - idum, m1);
    ix1 = fem::mod(ia1 * ix1 + ic1, m1);
    ix2 = fem::mod(ix1, m2);
    ix1 = fem::mod(ia1 * ix1 + ic1, m1);
    ix3 = fem::mod(ix1, m3);
    FEM_DO_SAFE(j, 1, 97) {
      ix1 = fem::mod(ia1 * ix1 + ic1, m1);
      ix2 = fem::mod(ia2 * ix2 + ic2, m2);
      r(j) = (fem::dble(ix1) + fem::dble(ix2) * rm2) * rm1;
    }
    idum = 1;
  }
  ix1 = fem::mod(ia1 * ix1 + ic1, m1);
  ix2 = fem::mod(ia2 * ix2 + ic2, m2);
  const int ia3 = 4561;
  const int ic3 = 51349;
  ix3 = fem::mod(ia3 * ix3 + ic3, m3);
  //Clin-7/01/02       j = 1 + (97 * i x 3) / m3
  j = 1 + (97 * ix3) / m3;
  //Clin-4/2008:
  //C      if (j .gt. 97 .or. j .lt. 1) pause
  if (j > 97 || j < 1) {
    write(6, star), "In zpc ran1, j<1 or j>97", j;
  }
  return_value = r(j);
  r(j) = (fem::dble(ix1) + fem::dble(ix2) * rm2) * rm1;
  //C
  //Clin-6/23/00 check random number generator:
  number++;
  //C      if(number.le.100000) write(99,*) 'number, ran1=', number,ran1
  //C
  return return_value;
}

struct posit1_save
{
  int iseed;

  posit1_save() :
    iseed(fem::int0)
  {}
};

void
posit1(
  common& cmn,
  double& x,
  double& y,
  double const& r0)
{
  FEM_CMN_SVE(posit1);
  int& iseed = sve.iseed;
  //C
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  statement_10:
  x = 2e0 * ran1(cmn, iseed) - 1e0;
  y = 2e0 * ran1(cmn, iseed) - 1e0;
  if (fem::pow2(x) + fem::pow2(y) > 1e0) {
    goto statement_10;
  }
  x = x * r0;
  y = y * r0;
  //C
}

struct posit2_save
{
  int iseed;

  posit2_save() :
    iseed(fem::int0)
  {}
};

void
posit2(
  common& cmn,
  double& x,
  double& y)
{
  FEM_CMN_SVE(posit2);
  // SAVE
  int& iseed = sve.iseed;
  //
  //C
  //Cc      SAVE /ilist3/
  //Cc      SAVE /rndm3/
  iseed = cmn.iseedp;
  x = 2e0 * ran1(cmn, iseed) - 1e0;
  y = 2e0 * ran1(cmn, iseed) - 1e0;
  x = x * 5e0 * cmn.size1;
  y = y * 5e0 * cmn.size2;
  //C
}

struct posit3_save
{
  int iseed;

  posit3_save() :
    iseed(fem::int0)
  {}
};

void
posit3(
  common& cmn,
  double& x,
  double& y,
  double& z)
{
  FEM_CMN_SVE(posit3);
  // SAVE
  int& iseed = sve.iseed;
  //
  //C
  //Cc      SAVE /ilist3/
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  x = 2e0 * ran1(cmn, iseed) - 1e0;
  y = 2e0 * ran1(cmn, iseed) - 1e0;
  z = 2e0 * ran1(cmn, iseed) - 1e0;
  x = x * 5e0 * cmn.size1;
  y = y * 5e0 * cmn.size2;
  z = z * 5e0 * cmn.size3;
  //C
}

struct energy_save
{
  int iseed;

  energy_save() :
    iseed(fem::int0)
  {}
};

void
energy(
  common& cmn,
  double& e,
  double const& temp)
{
  FEM_CMN_SVE(energy);
  int& iseed = sve.iseed;
  //C
  //C       to generate the magnitude of the momentum e,
  //C       knowing the temperature of the local thermal distribution temp
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  statement_1000:
  //C
  e = ran1(cmn, iseed);
  e = e * ran1(cmn, iseed);
  e = e * ran1(cmn, iseed);
  //C
  if (e <= 0e0) {
    goto statement_1000;
  }
  e = -temp * fem::log(e);
  if (ran1(cmn, iseed) > fem::exp((e - fem::dsqrt(fem::pow2(e) +
      fem::pow2(cmn.xmp))) / temp)) {
    goto statement_1000;
  }
  //C
}

struct momntm_save
{
  double cost;
  int iseed;
  double phi;
  double sint;

  momntm_save() :
    cost(fem::double0),
    iseed(fem::int0),
    phi(fem::double0),
    sint(fem::double0)
  {}
};

void
momntm(
  common& cmn,
  double& px,
  double& py,
  double& pz,
  double const& e)
{
  FEM_CMN_SVE(momntm);
  // SAVE
  double& cost = sve.cost;
  int& iseed = sve.iseed;
  double& phi = sve.phi;
  double& sint = sve.sint;
  //
  //C
  //C       to generate the 3 components of the momentum px, py, pz,
  //C       from the magnitude of the momentum e
  //C
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  cost = 2e0 * ran1(cmn, iseed) - 1e0;
  //C     7/20/01:
  //C        sint = sqrt(1d0 - cost ** 2)
  sint = fem::dsqrt(1e0 - fem::pow2(cost));
  const double pi = 3.14159265358979e0;
  phi = 2e0 * pi * ran1(cmn, iseed);
  //C
  px = e * sint * fem::cos(phi);
  py = e * sint * fem::sin(phi);
  pz = e * cost;
  //C
}

struct lorenz_save
{
  double beta2;
  double gam;

  lorenz_save() :
    beta2(fem::double0),
    gam(fem::double0)
  {}
};

//C
//C*****************************************************************************
//C
void
lorenz(
  common& cmn,
  double const& energy,
  double const& px,
  double const& py,
  double const& pz,
  double const& bex,
  double const& bey,
  double const& bez)
{
  FEM_CMN_SVE(lorenz);
  common_write write(cmn);
  // COMMON lor
  double& enenew = cmn.enenew;
  double& pxnew = cmn.pxnew;
  double& pynew = cmn.pynew;
  double& pznew = cmn.pznew;
  //
  // SAVE
  double& beta2 = sve.beta2;
  double& gam = sve.gam;
  //
  //C
  //C     add in a cut for beta2 to prevent gam to be nan (infinity)
  //C
  //Cc      SAVE /lor/
  //C
  beta2 = fem::pow2(bex) + fem::pow2(bey) + fem::pow2(bez);
  if (beta2 == 0e0) {
    enenew = energy;
    pxnew = px;
    pynew = py;
    pznew = pz;
  }
  else {
    if (beta2 > 0.999999999999999e0) {
      beta2 = 0.999999999999999e0;
      write(6, star), "beta2=0.999999999999999";
    }
    //Clin-7/20/01:
    //C         gam = 1.d0 / sqrt(1.d0 - beta2)
    gam = 1.e0 / fem::dsqrt(1.e0 - beta2);
    enenew = gam * (energy - bex * px - bey * py - bez * pz);
    pxnew = -gam * bex * energy + (1.e0 + (gam - 1.e0) * fem::pow2(
      bex) / beta2) * px + (gam - 1.e0) * bex * bey / beta2 * py + (
      gam - 1.e0) * bex * bez / beta2 * pz;
    pynew = -gam * bey * energy + (gam - 1.e0) * bex * bey / beta2 *
      px + (1.e0 + (gam - 1.e0) * fem::pow2(bey) / beta2) * py + (
      gam - 1.e0) * bey * bez / beta2 * pz;
    pznew = -gam * bez * energy + (gam - 1.e0) * bex * bez / beta2 *
      px + (gam - 1.e0) * bey * bez / beta2 * py + (1.e0 + (gam -
      1.e0) * fem::pow2(bez) / beta2) * pz;
  }
  //C
}

struct genei_save
{
  double bex;
  double bey;
  double bez;
  double deta;
  double e;
  double etamax;
  double etamin;
  int i;
  int incmul;
  int iseed;
  double px;
  double py;
  double pz;
  double r0;
  double tau0;
  double temp;
  double x;
  double y;
  double z;

  genei_save() :
    bex(fem::double0),
    bey(fem::double0),
    bez(fem::double0),
    deta(fem::double0),
    e(fem::double0),
    etamax(fem::double0),
    etamin(fem::double0),
    i(fem::int0),
    incmul(fem::int0),
    iseed(fem::int0),
    px(fem::double0),
    py(fem::double0),
    pz(fem::double0),
    r0(fem::double0),
    tau0(fem::double0),
    temp(fem::double0),
    x(fem::double0),
    y(fem::double0),
    z(fem::double0)
  {}
};

void
genei(
  common& cmn)
{
  FEM_CMN_SVE(genei);
  common_write write(cmn);
  // COMMON para1
  int& mul = cmn.mul;
  // COMMON para2
  double& xmp = cmn.xmp;
  // COMMON para5
  int& iconfg = cmn.iconfg;
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
  // COMMON prec5
  arr_ref<double> eta(cmn.eta, dimension(maxptn));
  //
  // SAVE
  double& bex = sve.bex;
  double& bey = sve.bey;
  double& bez = sve.bez;
  double& deta = sve.deta;
  double& e = sve.e;
  double& etamax = sve.etamax;
  double& etamin = sve.etamin;
  int& i = sve.i;
  int& incmul = sve.incmul;
  int& iseed = sve.iseed;
  double& px = sve.px;
  double& py = sve.py;
  double& pz = sve.pz;
  double& r0 = sve.r0;
  double& tau0 = sve.tau0;
  double& temp = sve.temp;
  double& x = sve.x;
  double& y = sve.y;
  double& z = sve.z;
  //
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para2/
  //Cc      SAVE /para3/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec1/
  //Cc      SAVE /prec5/
  //Cc      SAVE /lor/
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  incmul = 4000;
  temp = 0.5e0;
  etamin = -5e0;
  etamax = 5e0;
  r0 = 5e0;
  tau0 = 0.1e0;
  deta = etamax - etamin;
  //C
  FEM_DO_SAFE(i, mul + 1, mul + incmul) {
    ityp0(i) = 21;
    xmass0(i) = xmp;
    energy(cmn, e, temp);
    momntm(cmn, px, py, pz, e);
    //C     7/20/01:
    //C           e = sqrt(e ** 2 + xmp ** 2)
    e = fem::dsqrt(fem::pow2(e) + fem::pow2(xmp));
    if (iconfg <= 3) {
      eta(i) = etamin + deta * ran1(cmn, iseed);
      bex = 0e0;
      bey = 0e0;
      bez = -fem::tanh(eta(i));
      lorenz(cmn, e, px, py, pz, bex, bey, bez);
      px0(i) = cmn.pxnew;
      py0(i) = cmn.pynew;
      pz0(i) = cmn.pznew;
      e0(i) = cmn.enenew;
    }
    else {
      px0(i) = px;
      py0(i) = py;
      pz0(i) = pz;
      e0(i) = e;
    }
  }
  //C
  FEM_DO_SAFE(i, mul + 1, mul + incmul) {
    if (iconfg <= 3) {
      gz0(i) = tau0 * fem::sinh(eta(i));
      ft0(i) = tau0 * fem::cosh(eta(i));
      if (iconfg == 1) {
        posit1(cmn, x, y, r0);
        gx0(i) = x + px0(i) * ft0(i) / e0(i);
        gy0(i) = y + py0(i) * ft0(i) / e0(i);
      }
      else if (iconfg == 2 || iconfg == 3) {
        posit2(cmn, x, y);
        gx0(i) = x;
        gy0(i) = y;
      }
    }
    else {
      ft0(i) = 0e0;
      posit3(cmn, x, y, z);
      gx0(i) = x;
      gy0(i) = y;
      gz0(i) = z;
    }
  }
  //C
  mul += incmul;
  //C
  //C       check if it's necessary to adjust array size 'adarr'
  if (mul >= maxptn || mul == 0) {
    write(6, star), "event", cmn.ievt, "has", mul, "number of gluon",
      "adjusting counting is necessary";
    FEM_STOP("adarr");
  }
  //C
}

struct boosti_save
{
  double bex;
  double bey;
  double bez;
  double e1;
  int i;
  double px1;
  double py1;
  double pz1;

  boosti_save() :
    bex(fem::double0),
    bey(fem::double0),
    bez(fem::double0),
    e1(fem::double0),
    i(fem::int0),
    px1(fem::double0),
    py1(fem::double0),
    pz1(fem::double0)
  {}
};

void
boosti(
  common& cmn)
{
  FEM_CMN_SVE(boosti);
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
  // COMMON lor
  double& enenew = cmn.enenew;
  double& pxnew = cmn.pxnew;
  double& pynew = cmn.pynew;
  double& pznew = cmn.pznew;
  //
  // SAVE
  double& bex = sve.bex;
  double& bey = sve.bey;
  double& bez = sve.bez;
  double& e1 = sve.e1;
  int& i = sve.i;
  double& px1 = sve.px1;
  double& py1 = sve.py1;
  double& pz1 = sve.pz1;
  //
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para6/
  //Cc      SAVE /prec1/
  //Cc      SAVE /lor/
  //C
  bex = 0e0;
  bey = 0e0;
  bez = -fem::tanh(cmn.centy);
  //C
  //C       save data for many runs of the same initial condition
  FEM_DO_SAFE(i, 1, cmn.mul) {
    px1 = gx0(i);
    py1 = gy0(i);
    pz1 = gz0(i);
    e1 = ft0(i);
    lorenz(cmn, e1, px1, py1, pz1, bex, bey, bez);
    gx0(i) = pxnew;
    gy0(i) = pynew;
    gz0(i) = pznew;
    ft0(i) = enenew;
    px1 = px0(i);
    py1 = py0(i);
    pz1 = pz0(i);
    e1 = e0(i);
    lorenz(cmn, e1, px1, py1, pz1, bex, bey, bez);
    px0(i) = pxnew;
    py0(i) = pynew;
    pz0(i) = pznew;
    e0(i) = enenew;
  }
  //C
}

//C
//C*****************************************************************************
//C
void
inievt(
  common& cmn)
{
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para4/
  //C
  //Cbz1/25/99
  //C        mul = 0
  //Cbz1/25/99
  if (cmn.ireflg == 0) {
    readi(cmn);
  }
  if (cmn.igeflg != 0) {
    genei(cmn);
  }
  if (cmn.ibstfg != 0) {
    boosti(cmn);
  }
  //C
}

struct index1_save
{
  int i;
  int indxt;
  int ir;
  int j;
  int l;
  double q;

  index1_save() :
    i(fem::int0),
    indxt(fem::int0),
    ir(fem::int0),
    j(fem::int0),
    l(fem::int0),
    q(fem::double0)
  {}
};

void
index1(
  common& cmn,
  int const& n,
  int const& m,
  arr_cref<double> arrin,
  arr_ref<int> indx)
{
  FEM_CMN_SVE(index1);
  arrin(dimension(n));
  indx(dimension(n));
  int& i = sve.i;
  int& indxt = sve.indxt;
  int& ir = sve.ir;
  int& j = sve.j;
  int& l = sve.l;
  double& q = sve.q;
  //C     indexes the first m elements of ARRIN of length n, i.e., outputs INDX
  //C     such that ARRIN(INDEX(J)) is in ascending order for J=1,...,m
  //C
  FEM_DO_SAFE(j, 1, m) {
    indx(j) = j;
  }
  l = m / 2 + 1;
  ir = m;
  statement_10:
  if (l > 1) {
    l = l - 1;
    indxt = indx(l);
    q = arrin(indxt);
  }
  else {
    indxt = indx(ir);
    q = arrin(indxt);
    indx(ir) = indx(1);
    ir = ir - 1;
    if (ir == 1) {
      indx(1) = indxt;
      return;
    }
  }
  i = l;
  j = l + l;
  statement_20:
  if (j <= ir) {
    if (j < ir) {
      if (arrin(indx(j)) < arrin(indx(j + 1))) {
        j++;
      }
    }
    if (q < arrin(indx(j))) {
      indx(i) = indx(j);
      i = j;
      j += j;
    }
    else {
      j = ir + 1;
    }
    goto statement_20;
  }
  indx(i) = indxt;
  goto statement_10;
  //C
}

struct ftime1_save
{
  double aa;

  ftime1_save() :
    aa(fem::double0)
  {}
};

double
ftime1(
  common& cmn,
  int& iseed)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(ftime1);
  // SAVE
  double& aa = sve.aa;
  //
  //C
  //C       this program is used to generate formation time
  //C       the calling program needs a common /par1/
  //C       and declare external ftime1
  //C
  //Clin-8/19/02
  //C
  //Cc      SAVE /par1/
  //C
  const double hbarc = 0.197327054e0;
  aa = hbarc / cmn.formt;
  //C
  //Clin7/20/01:
  //C        ftime1 = aa * sqrt(1d0 / ran1(iseed) - 1d0)
  return_value = aa * fem::dsqrt(1e0 / ran1(cmn, iseed) - 1e0);
  return return_value;
}

struct ftime_save
{
  int i;
  int iseed;
  double xmt2;

  ftime_save() :
    i(fem::int0),
    iseed(fem::int0),
    xmt2(fem::double0)
  {}
};

void
ftime(
  common& cmn)
{
  FEM_CMN_SVE(ftime);
  int& mul = cmn.mul;
  const int maxptn = 400001;
  arr_ref<double> ft0(cmn.ft0, dimension(maxptn));
  arr_cref<double> px0(cmn.px0, dimension(maxptn));
  arr_cref<double> py0(cmn.py0, dimension(maxptn));
  arr_cref<double> e0(cmn.e0, dimension(maxptn));
  arr_ref<int> indx(cmn.indx, dimension(maxptn));
  arr_ref<double> ct(cmn.ct, dimension(maxptn));
  arr_ref<double> ot(cmn.ot, dimension(maxptn));
  double& tlarge = cmn.tlarge;
  int& isoft = cmn.isoft;
  //
  int& i = sve.i;
  int& iseed = sve.iseed;
  double& xmt2 = sve.xmt2;
  //C       this subroutine generates formation time for the particles
  //C       indexing ft(i)
  //C       input e(i)
  //C       output ft(i), indx(i)
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para2/
  //Cc      SAVE /para4/
  //Cc      SAVE /prec1/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //Cc      SAVE /par1/
  //Cc      SAVE /anim/
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  //Clin-6/07/02 initialize here to expedite compiling, instead in zpcbdt:
  FEM_DO_SAFE(i, 1, maxptn) {
    ct(i) = 0e0;
    ot(i) = 0e0;
  }
  tlarge = 1000000.e0;
  //Clin-6/07/02-end
  //C
  if (cmn.iftflg == 0) {
    //C     5/01/01 different prescription for parton initial formation time:
    if (isoft == 3 || isoft == 4 || isoft == 5) {
      FEM_DO_SAFE(i, 1, mul) {
        if (ft0(i) > tlarge) {
          ft0(i) = tlarge;
        }
      }
      goto statement_150;
    }
    else {
      //C     5/01/01-end
      //C
      FEM_DO_SAFE(i, 1, maxptn) {
        ft0(i) = tlarge;
      }
      FEM_DO_SAFE(i, 1, mul) {
        xmt2 = fem::pow2(px0(i)) + fem::pow2(py0(i)) + fem::pow2(cmn.xmp);
        cmn.formt = xmt2 / e0(i);
        ft0(i) = ftime1(cmn, iseed);
        if (ft0(i) > tlarge) {
          ft0(i) = tlarge;
        }
      }
      //C     5/01/01:
    }
    //C
  }
  //C
  //C     5/01/01:
  statement_150:
  //C
  //C        call index1(MAXPTN, mul, ft0, indx)
  if (mul > 1) {
    index1(cmn, maxptn, mul, ft0, indx);
  }
  else {
    //Clin-7/09/03: need to set value for mul=1:
    indx(1) = 1;
  }
  //C
}

struct inirec_save
{
  double energy;
  double formt;
  int i;
  int indxi;
  int iseed;

  inirec_save() :
    energy(fem::double0),
    formt(fem::double0),
    i(fem::int0),
    indxi(fem::int0),
    iseed(fem::int0)
  {}
};

void
inirec(
  common& cmn)
{
  FEM_CMN_SVE(inirec);
  common_write write(cmn);
  // COMMON para1
  int& mul = cmn.mul;
  // COMMON prec1
  const int maxptn = 400001;
  arr_cref<double> gx0(cmn.gx0, dimension(maxptn));
  arr_cref<double> gy0(cmn.gy0, dimension(maxptn));
  arr_cref<double> gz0(cmn.gz0, dimension(maxptn));
  arr_cref<double> ft0(cmn.ft0, dimension(maxptn));
  arr_cref<double> px0(cmn.px0, dimension(maxptn));
  arr_cref<double> py0(cmn.py0, dimension(maxptn));
  arr_cref<double> pz0(cmn.pz0, dimension(maxptn));
  arr_cref<double> e0(cmn.e0, dimension(maxptn));
  arr_cref<double> xmass0(cmn.xmass0, dimension(maxptn));
  arr_cref<int> ityp0(cmn.ityp0, dimension(maxptn));
  // COMMON prec2
  arr_ref<double> gx(cmn.gx, dimension(maxptn));
  arr_ref<double> gy(cmn.gy, dimension(maxptn));
  arr_ref<double> gz(cmn.gz, dimension(maxptn));
  arr_ref<double> ft(cmn.ft, dimension(maxptn));
  arr_ref<double> px(cmn.px, dimension(maxptn));
  arr_ref<double> py(cmn.py, dimension(maxptn));
  arr_ref<double> pz(cmn.pz, dimension(maxptn));
  arr_ref<double> e(cmn.e, dimension(maxptn));
  arr_ref<double> xmass(cmn.xmass, dimension(maxptn));
  arr_ref<int> ityp(cmn.ityp, dimension(maxptn));
  // COMMON prec3
  arr_ref<double> gxs(cmn.gxs, dimension(maxptn));
  arr_ref<double> gys(cmn.gys, dimension(maxptn));
  arr_ref<double> gzs(cmn.gzs, dimension(maxptn));
  arr_ref<double> fts(cmn.fts, dimension(maxptn));
  arr_ref<double> pxs(cmn.pxs, dimension(maxptn));
  arr_ref<double> pys(cmn.pys, dimension(maxptn));
  arr_ref<double> pzs(cmn.pzs, dimension(maxptn));
  arr_ref<double> es(cmn.es, dimension(maxptn));
  arr_ref<double> xmasss(cmn.xmasss, dimension(maxptn));
  arr_ref<int> ityps(cmn.ityps, dimension(maxptn));
  // COMMON prec4
  arr_ref<double> vx(cmn.vx, dimension(maxptn));
  arr_ref<double> vy(cmn.vy, dimension(maxptn));
  arr_ref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON prec5
  arr_ref<double> eta(cmn.eta, dimension(maxptn));
  arr_ref<double> rap(cmn.rap, dimension(maxptn));
  arr_ref<double> tau(cmn.tau, dimension(maxptn));
  // COMMON prec6
  arr_ref<double> etas(cmn.etas, dimension(maxptn));
  arr_ref<double> raps(cmn.raps, dimension(maxptn));
  arr_ref<double> taus(cmn.taus, dimension(maxptn));
  // COMMON ilist4
  arr_cref<int> indx(cmn.indx, dimension(maxptn));
  // COMMON ilist7
  arr_cref<int> lstrg0(cmn.lstrg0, dimension(maxptn));
  arr_cref<int> lpart0(cmn.lpart0, dimension(maxptn));
  // COMMON ilist8
  arr_ref<int> lstrg1(cmn.lstrg1, dimension(maxptn));
  arr_ref<int> lpart1(cmn.lpart1, dimension(maxptn));
  // COMMON precpb
  arr_ref<double> vxp(cmn.vxp, dimension(maxptn));
  arr_ref<double> vyp(cmn.vyp, dimension(maxptn));
  arr_ref<double> vzp(cmn.vzp, dimension(maxptn));
  // COMMON precpa
  arr_cref<double> vxp0(cmn.vxp0, dimension(maxptn));
  arr_cref<double> vyp0(cmn.vyp0, dimension(maxptn));
  arr_cref<double> vzp0(cmn.vzp0, dimension(maxptn));
  arr_ref<double> xstrg0(cmn.xstrg0, dimension(maxptn));
  arr_ref<double> ystrg0(cmn.ystrg0, dimension(maxptn));
  arr_cref<double> xstrg(cmn.xstrg, dimension(maxptn));
  arr_cref<double> ystrg(cmn.ystrg, dimension(maxptn));
  arr_ref<int> istrg0(cmn.istrg0, dimension(maxptn));
  arr_cref<int> istrg(cmn.istrg, dimension(maxptn));
  // COMMON anim
  int& isoft = cmn.isoft;
  // COMMON frzprc
  arr_ref<double> gxfrz(cmn.gxfrz, dimension(maxptn));
  arr_ref<double> gyfrz(cmn.gyfrz, dimension(maxptn));
  arr_ref<double> gzfrz(cmn.gzfrz, dimension(maxptn));
  arr_ref<double> ftfrz(cmn.ftfrz, dimension(maxptn));
  arr_ref<double> pxfrz(cmn.pxfrz, dimension(maxptn));
  arr_ref<double> pyfrz(cmn.pyfrz, dimension(maxptn));
  arr_ref<double> pzfrz(cmn.pzfrz, dimension(maxptn));
  arr_ref<double> efrz(cmn.efrz, dimension(maxptn));
  arr_ref<double> xmfrz(cmn.xmfrz, dimension(maxptn));
  arr_ref<int> ifrz(cmn.ifrz, dimension(maxptn));
  arr_ref<int> idfrz(cmn.idfrz, dimension(maxptn));
  // COMMON para7
  int& ioscar = cmn.ioscar;
  //
  // SAVE
  double& energy = sve.energy;
  double& formt = sve.formt;
  int& i = sve.i;
  int& indxi = sve.indxi;
  //
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para4/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec1/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec3/
  //Cc      SAVE /prec4/
  //Cc      SAVE /prec5/
  //Cc      SAVE /prec6/
  //Cc      SAVE /ilist4/
  //Cbz1/25/99
  //Cc      SAVE /ilist7/
  //Cc      SAVE /ilist8/
  //Cbz1/25/99end
  //Cc      SAVE /smearz/
  //Clin-8/2015:
  //C        dimension vxp(MAXPTN), vyp(MAXPTN), vzp(MAXPTN)
  //Clin-8/2015:
  //C        common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
  //Cc      SAVE /precpa/
  //Cc      SAVE /anim/
  //Clin-6/06/02 local parton freezeout:
  //Cc      SAVE /frzprc/
  //Cc      SAVE /rndm3/
  sve.iseed = cmn.iseedp;
  //Clin-6/06/02 local freezeout initialization:
  if (isoft == 5) {
    cmn.itlast = 0;
    inifrz();
  }
  //C
  FEM_DO_SAFE(i, 1, mul) {
    //Clin-7/09/01 define indx(i) to save time:
    //C           ityp(i) = ityp0(indx(i))
    //C           gx(i) = gx0(indx(i))
    //C           gy(i) = gy0(indx(i))
    //C           gz(i) = gz0(indx(i))
    //C           ft(i) = ft0(indx(i))
    //C           px(i) = px0(indx(i))
    //C           py(i) = py0(indx(i))
    //C           pz(i) = pz0(indx(i))
    //C           e(i) = e0(indx(i))
    //C           xmass(i) = xmass0(indx(i))
    //Ccbz1/25/99
    //C           LSTRG1(I) = LSTRG0(INDX(I))
    //C           LPART1(I) = LPART0(INDX(I))
    //Ccbz1/25/99end
    indxi = indx(i);
    ityp(i) = ityp0(indxi);
    gx(i) = gx0(indxi);
    gy(i) = gy0(indxi);
    gz(i) = gz0(indxi);
    ft(i) = ft0(indxi);
    px(i) = px0(indxi);
    py(i) = py0(indxi);
    pz(i) = pz0(indxi);
    e(i) = e0(indxi);
    xmass(i) = xmass0(indxi);
    lstrg1(i) = lstrg0(indxi);
    lpart1(i) = lpart0(indxi);
    vxp(i) = vxp0(indxi);
    vyp(i) = vyp0(indxi);
    vzp(i) = vzp0(indxi);
    //Clin-8/2015:
    xstrg0(i) = xstrg(indxi);
    ystrg0(i) = ystrg(indxi);
    istrg0(i) = istrg(indxi);
    //Clin-7/09/01-end
    //C
    //Clin-6/06/02 local freezeout initialization:
    if (isoft == 5) {
      idfrz(i) = ityp(i);
      gxfrz(i) = gx(i);
      gyfrz(i) = gy(i);
      gzfrz(i) = gz(i);
      ftfrz(i) = ft(i);
      pxfrz(i) = px(i);
      pyfrz(i) = py(i);
      pzfrz(i) = pz(i);
      efrz(i) = e(i);
      xmfrz(i) = xmass(i);
      ifrz(i) = 0;
    }
    //Clin-6/06/02-end
  }
  //C
  //C       save particle info for fixed time analysis
  FEM_DO_SAFE(i, 1, mul) {
    ityps(i) = ityp(i);
    gxs(i) = gx(i);
    gys(i) = gy(i);
    gzs(i) = gz(i);
    fts(i) = ft(i);
    pxs(i) = px(i);
    pys(i) = py(i);
    pzs(i) = pz(i);
    es(i) = e(i);
    xmasss(i) = xmass(i);
  }
  //C
  //Clin-6/2009
  if (isoft == 1 && (ioscar == 2 || ioscar == 3)) {
    write(92, star), cmn.iaevt, cmn.miss, mul;
  }
  //C
  FEM_DO_SAFE(i, 1, mul) {
    energy = e(i);
    vx(i) = px(i) / energy;
    vy(i) = py(i) / energy;
    vz(i) = pz(i) / energy;
    if (cmn.iftflg == 0) {
      formt = ft(i);
      //C     7/09/01 propagate partons with parent velocity till formation
      //C     so that partons in same hadron have 0 distance:
      //C            gx(i) = gx(i) + vx(i) * formt
      //C            gy(i) = gy(i) + vy(i) * formt
      //C            gz(i) = gz(i) + vz(i) * formt
      if (isoft == 3 || isoft == 4 || isoft == 5) {
        gx(i) += vxp(i) * formt;
        gy(i) += vyp(i) * formt;
        gz(i) += vzp(i) * formt;
      }
      else {
        gx(i) += vx(i) * formt;
        gy(i) += vy(i) * formt;
        gz(i) += vz(i) * formt;
      }
      //C     7/09/01-end
      //C
      //C     3/27/00-ctest off no smear z on partons to avoid eta overflow:
      //C              gz(i) = gz(i)+smearp*(2d0 * ran1(iseed) - 1d0)
      //C     to give eta=y +- smearp*random:
      //C              smeary=smearp*(2d0 * ran1(iseed) - 1d0)
      //C              smearf=dexp(2*smeary)*(1+vz(i))/(1-vz(i)+1.d-8)
      //C              gz(i) = gz(i)+formt*(smearf-1)/(smearf+1)
      //C     3/27/00-end
    }
    //C
    //Clin-6/2009 write out initial parton information after string melting
    //C     and after propagating to its format time:
    if (ioscar == 2 || ioscar == 3) {
      if (fem::dmax1(fem::abs(gx(i)), fem::abs(gy(i)), fem::abs(gz(i)),
          fem::abs(ft(i))) < 9999) {
        //Clin-8/2015:
        write(92,
          "(i3,2(1x,f7.2),1x,f8.2,1x,f6.3,4(1x,f8.2),1x,i5,2(1x,f7.2))"),
          ityp(i), px(i), py(i), pz(i), xmass(i), gx(i), gy(i), gz(i),
          ft(i), istrg0(i), xstrg0(i), ystrg0(i);
      }
      else {
        //Clin-8/2015:
        write(92,
          "(i3,2(1x,f7.2),1x,f8.2,1x,f6.3,4(1x,e8.2),1x,i5,2(1x,f7.2))"),
          ityp(i), px(i), py(i), pz(i), xmass(i), gx(i), gy(i), gz(i),
          ft(i), istrg0(i), xstrg0(i), ystrg0(i);
      }
    }
    //Clin-8/2015:
    //C 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
    //C 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
    //C     reduce file size:
    //C 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f9.3),
    //C     1          1x,I6,2(1x,f8.3))
    //C 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e9.3),
    //C     1          1x,I6,2(1x,f8.3))
    //C
  }
  //C
  if (cmn.iconfg <= 3) {
    FEM_DO_SAFE(i, 1, mul) {
      if (ft(i) <= fem::abs(gz(i))) {
        eta(i) = 1000000.e0;
      }
      else {
        eta(i) = 0.5e0 * fem::log((ft(i) + gz(i)) / (ft(i) - gz(i)));
      }
      if (e(i) <= fem::abs(pz(i))) {
        rap(i) = 1000000.e0;
      }
      else {
        rap(i) = 0.5e0 * fem::log((e(i) + pz(i)) / (e(i) - pz(i)));
      }
      //Clin-8/2015 to avoid IEEE_OVERFLOW_FLAG:
      //C              tau(i) = ft(i) / cosh(eta(i))
      if (eta(i) < 1000000.e0) {
        tau(i) = ft(i) / fem::cosh(eta(i));
      }
      else {
        tau(i) = 1e-10;
      }
      //C
    }
    //C
    FEM_DO_SAFE(i, 1, mul) {
      etas(i) = eta(i);
      raps(i) = rap(i);
      taus(i) = tau(i);
    }
  }
  //C
}

struct iilist_save
{
  int i;
  int i1;
  int i2;
  int i3;

  iilist_save() :
    i(fem::int0),
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0)
  {}
};

void
iilist(
  common& cmn)
{
  FEM_CMN_SVE(iilist);
  // COMMON para1
  int& mul = cmn.mul;
  // COMMON ilist1
  const int maxptn = 400001;
  arr_ref<int> next(cmn.next, dimension(maxptn));
  arr_ref<int> last(cmn.last, dimension(maxptn));
  arr_ref<int> icsta(cmn.icsta, dimension(maxptn));
  arr_ref<int> nic(cmn.nic, dimension(maxptn));
  arr_ref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist2
  arr_ref<int, 3> icel(cmn.icel, dimension(10, 10, 10));
  // COMMON ilist5
  arr_ref<double> ct(cmn.ct, dimension(maxptn));
  arr_ref<double> ot(cmn.ot, dimension(maxptn));
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  int& i = sve.i;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  //
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //Cc      SAVE /ilist6/
  //C
  cmn.iscat = maxptn;
  cmn.jscat = maxptn;
  //C
  FEM_DO_SAFE(i, 1, mul) {
    next(i) = 0;
    last(i) = 0;
    icsta(i) = 0;
    nic(i) = 0;
    icels(i) = 0;
  }
  //C
  cmn.icell = 0;
  FEM_DO_SAFE(i1, 1, 10) {
    FEM_DO_SAFE(i2, 1, 10) {
      FEM_DO_SAFE(i3, 1, 10) {
        icel(i1, i2, i3) = 0;
      }
    }
  }
  //C
  cmn.ichkpt = 0;
  cmn.ifmpt = 1;
  //C
  FEM_DO_SAFE(i, 1, mul) {
    ct(i) = tlarge;
    ot(i) = tlarge;
  }
  //C
  cmn.iopern = 0;
  cmn.icolln = 0;
  cmn.t = 0.e0;
  //C
}

struct inian2_save
{
  int i;

  inian2_save() :
    i(fem::int0)
  {}
};

void
inian2(
  common& cmn)
{
  FEM_CMN_SVE(inian2);
  // COMMON ana2
  arr_ref<double> det(cmn.det, dimension(12));
  arr_ref<double> dn(cmn.dn, dimension(12));
  arr_ref<double> det1(cmn.det1, dimension(12));
  arr_ref<double> dn1(cmn.dn1, dimension(12));
  arr_ref<double> det2(cmn.det2, dimension(12));
  arr_ref<double> dn2(cmn.dn2, dimension(12));
  //
  // SAVE
  int& i = sve.i;
  //
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /ana2/
  //C
  if (cmn.iconfg <= 3) {
    FEM_DO_SAFE(i, 1, 12) {
      det(i) = 0e0;
      dn(i) = 0e0;
      det1(i) = 0e0;
      dn1(i) = 0e0;
      det2(i) = 0e0;
      dn2(i) = 0e0;
    }
  }
  //C
}

//C
//C*****************************************************************************
//C
void
inirun(
  common& cmn)
{
  //C
  //C       sort prec2 according to increasing formation time
  ftime(cmn);
  inirec(cmn);
  iilist(cmn);
  inian2(cmn);
  //C
}

void
savrec(
  common& cmn,
  int const& i)
{
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  arr_cref<double> xmass(cmn.xmass, dimension(maxptn));
  arr_cref<int> ityp(cmn.ityp, dimension(maxptn));
  // COMMON prec3
  arr_ref<double> gxs(cmn.gxs, dimension(maxptn));
  arr_ref<double> gys(cmn.gys, dimension(maxptn));
  arr_ref<double> gzs(cmn.gzs, dimension(maxptn));
  arr_ref<double> fts(cmn.fts, dimension(maxptn));
  arr_ref<double> pxs(cmn.pxs, dimension(maxptn));
  arr_ref<double> pys(cmn.pys, dimension(maxptn));
  arr_ref<double> pzs(cmn.pzs, dimension(maxptn));
  arr_ref<double> es(cmn.es, dimension(maxptn));
  arr_ref<double> xmasss(cmn.xmasss, dimension(maxptn));
  arr_ref<int> ityps(cmn.ityps, dimension(maxptn));
  // COMMON prec5
  arr_cref<double> eta(cmn.eta, dimension(maxptn));
  arr_cref<double> rap(cmn.rap, dimension(maxptn));
  arr_cref<double> tau(cmn.tau, dimension(maxptn));
  // COMMON prec6
  arr_ref<double> etas(cmn.etas, dimension(maxptn));
  arr_ref<double> raps(cmn.raps, dimension(maxptn));
  arr_ref<double> taus(cmn.taus, dimension(maxptn));
  //
  //C
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec3/
  //Cc      SAVE /prec5/
  //Cc      SAVE /prec6/
  //C
  ityps(i) = ityp(i);
  gxs(i) = gx(i);
  gys(i) = gy(i);
  gzs(i) = gz(i);
  fts(i) = ft(i);
  pxs(i) = px(i);
  pys(i) = py(i);
  pzs(i) = pz(i);
  es(i) = e(i);
  xmasss(i) = xmass(i);
  etas(i) = eta(i);
  raps(i) = rap(i);
  taus(i) = tau(i);
  //C
}

struct getict_save
{
  int i;

  getict_save() :
    i(fem::int0)
  {}
};

void
getict(
  common& cmn,
  double& t1)
{
  FEM_CMN_SVE(getict);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON ilist1
  int& iscat = cmn.iscat;
  int& jscat = cmn.jscat;
  arr_cref<int> next(cmn.next, dimension(maxptn));
  int& ictype = cmn.ictype;
  arr_cref<int> icsta(cmn.icsta, dimension(maxptn));
  // COMMON ilist4
  int& ifmpt = cmn.ifmpt;
  // COMMON ilist5
  arr_cref<double> ot(cmn.ot, dimension(maxptn));
  //
  // SAVE
  int& i = sve.i;
  //
  //Cc      SAVE /para1/
  //Cc      SAVE /prec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //C
  //C       neglect possibility of 2 collisions at the same time
  //C0       set initial conditions
  //C
  t1 = cmn.tlarge;
  iscat = 0;
  jscat = 0;
  //C
  //C1      get next collision between particles
  FEM_DO_SAFE(i, 1, cmn.ichkpt) {
    if (ot(i) < t1) {
      t1 = ot(i);
      iscat = i;
    }
  }
  if (iscat != 0) {
    jscat = next(iscat);
  }
  //C
  //C2      get ictype
  //C     10/30/02 ictype=0:collision; 1:parton formation
  if (iscat != 0 && jscat != 0) {
    if (icsta(iscat) == 0 && icsta(jscat) == 0) {
      ictype = 0;
    }
    else {
      ictype = 4;
    }
  }
  else if (iscat != 0 || jscat != 0) {
    ictype = 3;
  }
  //C
  if (ifmpt <= cmn.mul) {
    if (ft(ifmpt) < t1) {
      ictype = 1;
      t1 = ft(ifmpt);
    }
    else if (ft(ifmpt) == t1) {
      if (ictype == 0) {
        ictype = 2;
      }
      if (ictype == 3) {
        ictype = 5;
      }
      if (ictype == 4) {
        ictype = 6;
      }
    }
  }
  //C
}

int
integ(
  double const& x)
{
  int return_value = fem::int0;
  //C       this function is used to get the largest integer that is smaller than
  //C       x
  //C
  if (x < 0e0) {
    return_value = fem::fint(x - 1e0);
  }
  else {
    return_value = fem::fint(x);
  }
  //C
  return return_value;
}

struct newcre_save
{
  int j;

  newcre_save() :
    j(fem::int0)
  {}
};

void
newcre(
  common& cmn,
  int const& i,
  int& k)
{
  FEM_CMN_SVE(newcre);
  // COMMON ilist1
  const int maxptn = 400001;
  arr_ref<int> nic(cmn.nic, dimension(maxptn));
  //
  // SAVE
  int& j = sve.j;
  //
  //C       this subroutine is used to mk rearrange of the new cell a particle
  //C       enters,
  //C       input i
  //C       output nic(i)
  //C
  //Cc      SAVE /ilist1/
  //C
  if (k == 0) {
    k = i;
    nic(i) = 0;
  }
  else if (nic(k) == 0) {
    nic(k) = i;
    nic(i) = k;
  }
  else {
    j = k;
    while (nic(j) != k) {
      j = nic(j);
    }
    //C
    nic(j) = i;
    nic(i) = k;
    //C
  }
  //C
}

struct celasn_save
{
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  double td;
  double tt;

  celasn_save() :
    i(fem::int0),
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    j(fem::int0),
    td(fem::double0),
    tt(fem::double0)
  {}
};

void
celasn(
  common& cmn)
{
  FEM_CMN_SVE(celasn);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_ref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist2
  int& icell = cmn.icell;
  arr_ref<int, 3> icel(cmn.icel, dimension(10, 10, 10));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  double& v1 = cmn.v1;
  double& v2 = cmn.v2;
  double& v3 = cmn.v3;
  //
  // SAVE
  int& i = sve.i;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& j = sve.j;
  double& td = sve.td;
  double& tt = sve.tt;
  //
  //C       this subroutine is used to assign a cell for a newly formed particle
  //C       output: nic(MAXPTN) icels(MAXPTN) in the common /ilist1/
  //C       icell, and icel(10,10,10) in the common /ilist2/
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist4/
  //C
  i = cmn.ifmpt;
  tt = ft(i);
  td = tt - cmn.size;
  if (iconfg == 1 && (size1 == 0e0 || size2 == 0e0 || size3 == 0e0)) {
    i1 = 11;
    i2 = 11;
    i3 = 11;
  }
  else if (iconfg == 4 || td <= 0e0) {
    i1 = integ(gx(i) / size1) + 6;
    i2 = integ(gy(i) / size2) + 6;
    i3 = integ(gz(i) / size3) + 6;
    if (integ(gx(i) / size1) == gx(i) / size1 && vx(i) < 0e0) {
      i1 = i1 - 1;
    }
    if (integ(gy(i) / size2) == gy(i) / size2 && vy(i) < 0e0) {
      i2 = i2 - 1;
    }
    if (integ(gz(i) / size3) == gz(i) / size3 && vz(i) < 0e0) {
      i3 = i3 - 1;
    }
  }
  else {
    i1 = integ(gx(i) / (size1 + v1 * td)) + 6;
    i2 = integ(gy(i) / (size2 + v2 * td)) + 6;
    i3 = integ(gz(i) / (size3 + v3 * td)) + 6;
    if (integ(gx(i) / (size1 + v1 * td)) == gx(i) / (size1 + v1 *
        td) && vx(i) < (i1 - 6) * v1) {
      i1 = i1 - 1;
    }
    if (integ(gy(i) / (size2 + v2 * td)) == gy(i) / (size2 + v2 *
        td) && vy(i) < (i2 - 6) * v2) {
      i2 = i2 - 1;
    }
    if (integ(gz(i) / (size3 + v3 * td)) == gz(i) / (size3 + v3 *
        td) && vz(i) < (i3 - 6) * v3) {
      i3 = i3 - 1;
    }
  }
  //C
  if (i1 <= 0 || i1 >= 11 || i2 <= 0 || i2 >= 11 || i3 <= 0 || i3 >= 11) {
    i1 = 11;
    i2 = 11;
    i3 = 11;
  }
  //C
  if (i1 == 11) {
    j = icell;
    newcre(cmn, i, j);
    icell = j;
    icels(i) = 111111;
  }
  else {
    j = icel(i1, i2, i3);
    newcre(cmn, i, j);
    icel(i1, i2, i3) = j;
    icels(i) = i1 * 10000 + i2 * 100 + i3;
  }
  //C
}

struct oldcre_save
{
  int j;

  oldcre_save() :
    j(fem::int0)
  {}
};

void
oldcre(
  common& cmn,
  int const& i)
{
  FEM_CMN_SVE(oldcre);
  // COMMON ilist1
  const int maxptn = 400001;
  arr_ref<int> nic(cmn.nic, dimension(maxptn));
  //
  // SAVE
  int& j = sve.j;
  //
  //C       this subroutine is used to rearrange the old cell nic when a particle
  //C       goes out of the cell
  //C
  //Cc      SAVE /ilist1/
  //C
  if (nic(i) == 0) {
    return;
  }
  //C
  j = nic(i);
  //C
  if (nic(j) == i) {
    nic(j) = 0;
    return;
  }
  //C
  while (nic(j) != i) {
    j = nic(j);
  }
  //C
  nic(j) = nic(i);
  //C
}

struct wallc1_save
{
  double t1;
  double t2;
  double t3;
  double tf;
  double v1p;
  double v2p;
  double v3p;
  double x1p;
  double x2p;
  double x3p;

  wallc1_save() :
    t1(fem::double0),
    t2(fem::double0),
    t3(fem::double0),
    tf(fem::double0),
    v1p(fem::double0),
    v2p(fem::double0),
    v3p(fem::double0),
    x1p(fem::double0),
    x2p(fem::double0),
    x3p(fem::double0)
  {}
};

void
wallc1(
  common& cmn,
  int const& i,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t,
  double& tmin)
{
  FEM_CMN_SVE(wallc1);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_ref<int> icsta(cmn.icsta, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  double& v1 = cmn.v1;
  double& v2 = cmn.v2;
  double& v3 = cmn.v3;
  double& size = cmn.size;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  double& t3 = sve.t3;
  double& tf = sve.tf;
  double& v1p = sve.v1p;
  double& v2p = sve.v2p;
  double& v3p = sve.v3p;
  double& x1p = sve.x1p;
  double& x2p = sve.x2p;
  double& x3p = sve.x3p;
  //
  //C       this subroutine is used to get wall collision time
  //C       when particle is inside the cube, it sets the icsta at the same time
  //C       input i,i1,i2,i3,t
  //C       output tmin, icsta(i)
  //C       note the icsta is not finally set. we need further judgement in
  //C       fixtim
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  x1p = gx(i);
  x2p = gy(i);
  x3p = gz(i);
  tf = ft(i);
  v1p = vx(i);
  v2p = vy(i);
  v3p = vz(i);
  //C
  if (t < size && tf < size) {
    //C
    if (v1p > 0e0) {
      t1 = ((fem::dble(i1) - 5e0) * size1 - x1p) / v1p + tf;
    }
    else if (v1p < 0e0) {
      t1 = ((fem::dble(i1) - 6e0) * size1 - x1p) / v1p + tf;
    }
    else {
      t1 = tlarge;
    }
    //C
    if (v2p > 0e0) {
      t2 = ((fem::dble(i2) - 5e0) * size2 - x2p) / v2p + tf;
    }
    else if (v2p < 0e0) {
      t2 = ((fem::dble(i2) - 6e0) * size2 - x2p) / v2p + tf;
    }
    else {
      t2 = tlarge;
    }
    //C
    if (v3p > 0e0) {
      t3 = ((fem::dble(i3) - 5e0) * size3 - x3p) / v3p + tf;
    }
    else if (v3p < 0e0) {
      t3 = ((fem::dble(i3) - 6e0) * size3 - x3p) / v3p + tf;
    }
    else {
      t3 = tlarge;
    }
    //C
    //C       if a particle is on the wall, we don't collide it on the same wall
    //C
    //C        if (t1 .eq. 0d0) t1 = tlarge
    //C        if (t2 .eq. 0d0) t2 = tlarge
    //C        if (t3 .eq. 0d0) t3 = tlarge
    //C
    tmin = fem::min(t1, t2, t3);
    //C
    //C       set icsta,
    //C       after checking this is not an earlier collision comparing with
    //C       a collision with another particle, we need to set icsta=0
    //C       after checking whether there is also a particle collision
    //C       at the same time, we need to reset the second bit of icsta
    //C
    if (tmin == t1) {
      if (v1p > 0e0) {
        icsta(i) = 101;
      }
      else {
        icsta(i) = 102;
      }
    }
    //C
    if (tmin == t2) {
      if (v2p > 0e0) {
        icsta(i) = 103;
      }
      else {
        icsta(i) = 104;
      }
    }
    //C
    if (tmin == t3) {
      if (v3p > 0e0) {
        icsta(i) = 105;
      }
      else {
        icsta(i) = 106;
      }
    }
    //C
    if (tmin <= size) {
      return;
    }
    //C
  }
  //C
  if (v1p > (i1 - 5) * v1) {
    t1 = ((i1 - 5) * (size1 - v1 * size) + v1p * tf - x1p) / (v1p - (
      i1 - 5) * v1);
  }
  else if (v1p < (i1 - 6) * v1) {
    t1 = ((i1 - 6) * (size1 - v1 * size) + v1p * tf - x1p) / (v1p - (
      i1 - 6) * v1);
  }
  else {
    t1 = tlarge;
  }
  //C
  if (v2p > (i2 - 5) * v2) {
    t2 = ((i2 - 5) * (size2 - v2 * size) + v2p * tf - x2p) / (v2p - (
      i2 - 5) * v2);
  }
  else if (v2p < (i2 - 6) * v2) {
    t2 = ((i2 - 6) * (size2 - v2 * size) + v2p * tf - x2p) / (v2p - (
      i2 - 6) * v2);
  }
  else {
    t2 = tlarge;
  }
  //C
  if (v3p > (i3 - 5) * v3) {
    t3 = ((i3 - 5) * (size3 - v3 * size) + v3p * tf - x3p) / (v3p - (
      i3 - 5) * v3);
  }
  else if (v3p < (i3 - 6) * v3) {
    t3 = ((i3 - 6) * (size3 - v3 * size) + v3p * tf - x3p) / (v3p - (
      i3 - 6) * v3);
  }
  else {
    t3 = tlarge;
  }
  //C
  //C       if a particle is on the wall, we don't collide it on the same wall
  //C
  //C        if (t1 .eq. 0d0) t1 = tlarge
  //C        if (t2 .eq. 0d0) t2 = tlarge
  //C        if (t3 .eq. 0d0) t3 = tlarge
  //C
  tmin = fem::min(t1, t2, t3);
  //C
  //C       set icsta,
  //C       after checking this is not an earlier collision comparing with
  //C       a collision with another particle, we need to set icsta=0
  //C       after checking whether there is also a particle collision
  //C       at the same time, we need to reset the second bit of icsta
  //C
  if (tmin == t1) {
    if (v1p > (i1 - 5) * v1) {
      icsta(i) = 101;
    }
    else {
      icsta(i) = 102;
    }
  }
  //C
  if (tmin == t2) {
    if (v2p > (i2 - 5) * v2) {
      icsta(i) = 103;
    }
    else {
      icsta(i) = 104;
    }
  }
  //C
  if (tmin == t3) {
    if (v3p > (i3 - 5) * v3) {
      icsta(i) = 105;
    }
    else {
      icsta(i) = 106;
    }
  }
  //C
}

struct wallc2_save
{
  double t1;
  double t2;
  double t3;
  double tf;
  double v1p;
  double v2p;
  double v3p;
  double x1p;
  double x2p;
  double x3p;

  wallc2_save() :
    t1(fem::double0),
    t2(fem::double0),
    t3(fem::double0),
    tf(fem::double0),
    v1p(fem::double0),
    v2p(fem::double0),
    v3p(fem::double0),
    x1p(fem::double0),
    x2p(fem::double0),
    x3p(fem::double0)
  {}
};

void
wallc2(
  common& cmn,
  int const& i,
  int const& /* i1 */,
  int const& /* i2 */,
  int const& /* i3 */,
  double const& /* t */,
  double& tmin)
{
  FEM_CMN_SVE(wallc2);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_ref<int> icsta(cmn.icsta, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  double& t3 = sve.t3;
  double& tf = sve.tf;
  double& v1p = sve.v1p;
  double& v2p = sve.v2p;
  double& v3p = sve.v3p;
  double& x1p = sve.x1p;
  double& x2p = sve.x2p;
  double& x3p = sve.x3p;
  //
  //C       this subroutine is used to get wall collision time
  //C       when particle is inside the cube, it sets the icsta at the same time
  //C       input i,i1,i2,i3,t
  //C       output tmin, icsta(i)
  //C       note the icsta is not finally set. we need further judgement in
  //C       fixtim
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  x1p = gx(i);
  x2p = gy(i);
  x3p = gz(i);
  tf = ft(i);
  v1p = vx(i);
  v2p = vy(i);
  v3p = vz(i);
  //C
  if (v1p > 0e0) {
    t1 = (5e0 * size1 - x1p) / v1p + tf;
  }
  else if (v1p < 0e0) {
    t1 = (-5e0 * size1 - x1p) / v1p + tf;
  }
  else {
    t1 = tlarge;
  }
  //C
  if (v2p > 0e0) {
    t2 = (5e0 * size2 - x2p) / v2p + tf;
  }
  else if (v2p < 0e0) {
    t2 = (-5e0 * size2 - x2p) / v2p + tf;
  }
  else {
    t2 = tlarge;
  }
  //C
  if (cmn.iconfg == 5) {
    if (v3p > 0e0) {
      t3 = (5e0 * size3 - x3p) / v3p + tf;
    }
    else if (v3p < 0e0) {
      t3 = (-5e0 * size3 - x3p) / v3p + tf;
    }
    else {
      t3 = tlarge;
    }
  }
  else {
    t3 = tlarge;
  }
  //C
  tmin = fem::min(t1, t2, t3);
  //C
  //C       set icsta,
  //C       after checking this is not an earlier collision comparing with
  //C       a collision with another particle, we need to set icsta=0
  //C       after checking whether there is also a particle collision
  //C       at the same time, we need to reset the second bit of icsta
  //C
  if (tmin == t1) {
    if (v1p > 0e0) {
      icsta(i) = 101;
    }
    else {
      icsta(i) = 102;
    }
  }
  //C
  if (tmin == t2) {
    if (v2p > 0e0) {
      icsta(i) = 103;
    }
    else {
      icsta(i) = 104;
    }
  }
  //C
  if (tmin == t3) {
    if (v3p > 0e0) {
      icsta(i) = 105;
    }
    else {
      icsta(i) = 106;
    }
  }
  //C
}

struct wallcb_save
{
  int icsta1;
  int icsta2;
  int icsta3;
  double t1;
  double t2;
  double t3;
  double tf;
  double v1p;
  double v2p;
  double v3p;
  double x1p;
  double x1pp;
  double x1q;
  double x2p;
  double x2pp;
  double x2q;
  double x3p;
  double x3pp;
  double x3q;

  wallcb_save() :
    icsta1(fem::int0),
    icsta2(fem::int0),
    icsta3(fem::int0),
    t1(fem::double0),
    t2(fem::double0),
    t3(fem::double0),
    tf(fem::double0),
    v1p(fem::double0),
    v2p(fem::double0),
    v3p(fem::double0),
    x1p(fem::double0),
    x1pp(fem::double0),
    x1q(fem::double0),
    x2p(fem::double0),
    x2pp(fem::double0),
    x2q(fem::double0),
    x3p(fem::double0),
    x3pp(fem::double0),
    x3q(fem::double0)
  {}
};

void
wallcb(
  common& cmn,
  int const& i,
  double const& t,
  double& tmin)
{
  FEM_CMN_SVE(wallcb);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_ref<int> icsta(cmn.icsta, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  double& v1 = cmn.v1;
  double& v2 = cmn.v2;
  double& v3 = cmn.v3;
  double& size = cmn.size;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  int& icsta1 = sve.icsta1;
  int& icsta2 = sve.icsta2;
  int& icsta3 = sve.icsta3;
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  double& t3 = sve.t3;
  double& tf = sve.tf;
  double& v1p = sve.v1p;
  double& v2p = sve.v2p;
  double& v3p = sve.v3p;
  double& x1p = sve.x1p;
  double& x1pp = sve.x1pp;
  double& x1q = sve.x1q;
  double& x2p = sve.x2p;
  double& x2pp = sve.x2pp;
  double& x2q = sve.x2q;
  double& x3p = sve.x3p;
  double& x3pp = sve.x3pp;
  double& x3q = sve.x3q;
  //
  //C       this subroutine is used to calculate the wall collision time
  //C       when the particle is outside the cube
  //C       input i,t
  //C       output tmin,icsta(i)
  //C       note the icsta is not finally set. we need further judgement in
  //C       fixtim
  //C
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       check if there is a collision by looking at the closest approach point
  //C       and see if it's inside the cube
  //C
  if (size1 == 0e0 || size2 == 0e0 || size3 == 0e0) {
    return;
  }
  //C
  x1p = gx(i);
  x2p = gy(i);
  x3p = gz(i);
  v1p = vx(i);
  v2p = vy(i);
  v3p = vz(i);
  tf = ft(i);
  //C
  if (t < size && tf < size) {
    if (x1p <  - 5e0 * size1 && v1p > 0e0) {
      t1 = (-5e0 * size1 - x1p) / v1p + tf;
    }
    else if (x1p > 5e0 * size1 && v1p < 0e0) {
      t1 = -(x1p - 5e0 * size1) / v1p + tf;
    }
    else {
      t1 = tlarge;
    }
    //C
    if (t1 != tlarge) {
      x2pp = x2p + v2p * (t1 - tf);
      x3pp = x3p + v3p * (t1 - tf);
      if (x2pp <=  - 5e0 * size2 || x2pp >= 5e0 * size2 || x3pp <=  -
          5e0 * size3 || x3pp >= 5e0 * size3) {
        t1 = tlarge;
      }
    }
    //C
    if (x2p <  - 5e0 * size2 && v2p > 0e0) {
      t2 = (-5e0 * size2 - x2p) / v2p + tf;
    }
    else if (x2p > 5e0 * size2 && v2p < 0e0) {
      t2 = -(x2p - 5e0 * size2) / v2p + tf;
    }
    else {
      t2 = tlarge;
    }
    //C
    if (t2 != tlarge) {
      x1pp = x1p + v1p * (t2 - tf);
      x3pp = x3p + v3p * (t2 - tf);
      if (x1pp <=  - 5e0 * size1 || x1pp >= 5e0 * size1 || x3pp <=  -
          5e0 * size3 || x3pp >= 5e0 * size3) {
        t2 = tlarge;
      }
    }
    //C
    if (x3p <  - 5e0 * size3 && v3p > 0e0) {
      t3 = (-5e0 * size3 - x3p) / v3p + tf;
    }
    else if (x3p > 5e0 * size3 && v3p < 0e0) {
      t3 = -(x3p - 5e0 * size3) / v3p + tf;
    }
    else {
      t3 = tlarge;
    }
    //C
    if (t3 != tlarge) {
      x1pp = x1p + v1p * (t3 - tf);
      x2pp = x2p + v2p * (t3 - tf);
      if (x1pp <=  - 5e0 * size1 || x1pp >= 5e0 * size1 || x2pp <=  -
          5e0 * size2 || x2pp >= 5e0 * size2) {
        t3 = tlarge;
      }
    }
    //C
    tmin = fem::min(t1, t2, t3);
    //C
    //C       set icsta,
    //C       after checking this is not an earlier collision comparing with
    //C       a collision with another particle, we need to set icsta=0
    //C       after checking whether there is also a particle collision
    //C       at the same time, we need to reset the second bit of icsta
    //C
    if (tmin == t1) {
      if (v1p > 0e0) {
        icsta(i) = 101;
      }
      else {
        icsta(i) = 102;
      }
    }
    //C
    if (tmin == t2) {
      if (v2p > 0e0) {
        icsta(i) = 103;
      }
      else {
        icsta(i) = 104;
      }
    }
    //C
    if (tmin == t3) {
      if (v3p > 0e0) {
        icsta(i) = 105;
      }
      else {
        icsta(i) = 106;
      }
    }
    //C
    if (tmin <= size) {
      return;
    }
    //C
  }
  //C
  //C       notice now x1q, x2q, x3q are coordinates at time t
  x1q = x1p + v1p * (t - tf);
  x2q = x2p + v2p * (t - tf);
  x3q = x3p + v3p * (t - tf);
  //C
  if (x1q <  - 5e0 * (size1 + v1 * (t - size)) && v1p >  - 5e0 * v1) {
    t1 = (-5e0 * (size1 - v1 * size) + v1p * tf - x1p) / (v1p - (-5e0) * v1);
    icsta1 = 101;
  }
  else if (x1q > 5e0 * (size1 + v1 * (t - size)) && v1p < 5e0 * v1) {
    t1 = (5e0 * (size1 - v1 * size) + v1p * tf - x1p) / (v1p - 5e0 * v1);
    icsta1 = 102;
  }
  else {
    t1 = tlarge;
  }
  //C
  if (t1 != tlarge) {
    x2pp = x2p + v2p * (t1 - tf);
    x3pp = x3p + v3p * (t1 - tf);
    if (x2pp <=  - 5e0 * (size2 + v2 * (t1 - size)) || x2pp >= 5e0 * (
        size2 + v2 * (t1 - size)) || x3pp <=  - 5e0 * (size3 + v3 * (
        t1 - size)) || x3pp >= 5e0 * (size3 + v3 * (t1 - size))) {
      t1 = tlarge;
    }
  }
  //C
  if (x2q <  - 5e0 * (size2 + v2 * (t - size)) && v2p >  - 5e0 * v2) {
    t2 = (-5e0 * (size2 - v2 * size) + v2p * tf - x2p) / (v2p - (-5e0) * v2);
    icsta2 = 103;
  }
  else if (x2q > 5e0 * (size2 + v2 * (t - size)) && v2p < 5e0 * v2) {
    t2 = (5e0 * (size2 - v2 * size) + v2p * tf - x2p) / (v2p - 5e0 * v2);
    icsta2 = 104;
  }
  else {
    t2 = tlarge;
  }
  //C
  if (t2 != tlarge) {
    x1pp = x1p + v1p * (t2 - tf);
    x3pp = x3p + v3p * (t2 - tf);
    if (x1pp <=  - 5e0 * (size1 + v1 * (t2 - size)) || x1pp >= 5e0 * (
        size1 + v1 * (t2 - size)) || x3pp <=  - 5e0 * (size3 + v3 * (
        t2 - size)) || x3pp >= 5e0 * (size3 + v3 * (t2 - size))) {
      t2 = tlarge;
    }
  }
  //C
  if (x3q <  - 5e0 * (size3 + v3 * (t - size)) && v3p >  - 5e0 * v3) {
    t3 = (-5e0 * (size3 - v3 * size) + v3p * tf - x3p) / (v3p - (-5e0) * v3);
    icsta3 = 105;
  }
  else if (x3q > 5e0 * (size3 + v3 * (t - size)) && v3p < 5e0 * v3) {
    t3 = (5e0 * (size3 - v3 * size) + v3p * tf - x3p) / (v3p - 5e0 * v3);
    icsta3 = 106;
  }
  else {
    t3 = tlarge;
  }
  //C
  if (t3 != tlarge) {
    x2pp = x2p + v2p * (t3 - tf);
    x1pp = x1p + v1p * (t3 - tf);
    if (x2pp <=  - 5e0 * (size2 + v2 * (t3 - size)) || x2pp >= 5e0 * (
        size2 + v2 * (t3 - size)) || x1pp <=  - 5e0 * (size1 + v1 * (
        t3 - size)) || x1pp >= 5e0 * (size1 + v1 * (t3 - size))) {
      t3 = tlarge;
    }
  }
  //C
  tmin = fem::min(t1, t2, t3);
  //C
  //C       set icsta,
  //C       after checking this is not an earlier collision comparing with
  //C       a collision with another particle, we need to set icsta=0
  //C       after checking whether there is also a particle collision
  //C       at the same time, we need to reset the second bit of icsta
  //C
  if (tmin == t1) {
    icsta(i) = icsta1;
  }
  else if (tmin == t2) {
    icsta(i) = icsta2;
  }
  else if (tmin == t3) {
    icsta(i) = icsta3;
  }
  //C
}

void
wallc(
  common& cmn,
  int const& i,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t,
  double& tmin)
{
  // COMMON para5
  int& iconfg = cmn.iconfg;
  //
  //C       this subroutine calculates the next time for collision with wall
  //C       for particle i
  //C       input particle label i,t
  //C       output tmin collision time with wall, icsta(i) wall collision
  //C       information
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /ilist5/
  //C
  tmin = cmn.tlarge;
  //C
  if (iconfg <= 2 || iconfg == 4) {
    //C       if particle is inside the cube
    if ((i1 >= 1 && i1 <= 10) || (i2 >= 1 && i2 <= 10) || (i3 >= 1 &&
        i3 <= 10)) {
      wallc1(cmn, i, i1, i2, i3, t, tmin);
      //C       if particle is outside the cube
    }
    else {
      wallcb(cmn, i, t, tmin);
    }
  }
  else if (iconfg == 3 || iconfg == 5) {
    wallc2(cmn, i, i1, i2, i3, t, tmin);
  }
  //C
}

struct fixtim_save
{
  int k;

  fixtim_save() :
    k(fem::int0)
  {}
};

void
fixtim(
  common& cmn,
  int const& l,
  double const& /* t */,
  double const& tmin1,
  double const& tmin,
  int const& nc)
{
  FEM_CMN_SVE(fixtim);
  // COMMON ilist1
  const int maxptn = 400001;
  arr_ref<int> next(cmn.next, dimension(maxptn));
  arr_ref<int> icsta(cmn.icsta, dimension(maxptn));
  // COMMON ilist5
  arr_cref<double> ct(cmn.ct, dimension(maxptn));
  arr_ref<double> ot(cmn.ot, dimension(maxptn));
  //
  // SAVE
  int& k = sve.k;
  //
  //C       this subroutine is used to compare the collision time with wall tmin1
  //C       and new collision time with particles for particle l
  //C       when used in ulist, input nc may be 0, which indicates no particle
  //C       collisions happen before wall collision, of course, then tmin=tmin1
  //C
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist5/
  //C
  k = nc;
  if (tmin < tmin1) {
    ot(l) = tmin;
    if (ct(l) < tmin1) {
      icsta(l) = 0;
    }
    else {
      icsta(l) += 10;
    }
    next(l) = k;
  }
  else if (tmin == tmin1) {
    ot(l) = tmin;
    if (nc == 0) {
      next(l) = 0;
    }
    else {
      icsta(l) += 10;
      next(l) = k;
    }
  }
  else {
    ot(l) = tmin1;
    next(l) = 0;
  }
  //C
}

struct isco1_save
{
  double a;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  int i1;
  int i2;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc2;
  double vp;

  isco1_save() :
    a(fem::double0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc2(fem::double0),
    vp(fem::double0)
  {}
};

void
isco1(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco1);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc2 = sve.tc2;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  i1 = i;
  i2 = j;
  //C
  p4 = ft(i2) - ft(i1);
  p1 = gx(i2) - gx(i1);
  p2 = gy(i2) - gy(i1);
  p3 = gz(i2) - gz(i1);
  //C
  q4 = e(i1);
  q1 = px(i1);
  q2 = py(i1);
  q3 = pz(i1);
  //C
  r4 = e(i2);
  r1 = px(i2);
  r2 = py(i2);
  r3 = pz(i2);
  //C
  a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
  b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
  c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
  d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
  ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
  f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
  //C
  //C       make sure particle 2 formed early
  h = a + b;
  if (h > 0e0) {
    g = a;
    a = -b;
    b = -g;
    //C
    g = c;
    c = d;
    d = g;
    //C
    i1 = j;
    i2 = i;
  }
  //C
  //C       check the approaching criteria
  if (allok) {
    //C
    vp = a * d - b * ee;
    //C
    allok = allok && vp < 0e0;
    //C
  }
  //C
  //C       check the closest approach distance criteria
  if (allok) {
    //C
    dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 * a * b *
      ee) / (fem::pow2(ee) - c * d);
    //C
    allok = allok && dm2 < cmn.cutof2;
    //C
  }
  //C
  //C       check the time criteria
  if (allok) {
    //C
    tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
    tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
    tm = 0.5e0 * (tc1 + tc2);
    //C
    allok = allok && tm > ft(i) && tm > ft(j);
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (h > 0e0) {
    t1 = tm;
    t2 = tm;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco2_save
{
  double a;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  int i1;
  int i2;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc2;
  double vp;

  isco2_save() :
    a(fem::double0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc2(fem::double0),
    vp(fem::double0)
  {}
};

void
isco2(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco2);
  // COMMON para5
  int& iordsc = cmn.iordsc;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc2 = sve.tc2;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  i1 = i;
  i2 = j;
  //C
  p4 = ft(i2) - ft(i1);
  p1 = gx(i2) - gx(i1);
  p2 = gy(i2) - gy(i1);
  p3 = gz(i2) - gz(i1);
  //C
  q4 = e(i1);
  q1 = px(i1);
  q2 = py(i1);
  q3 = pz(i1);
  //C
  r4 = e(i2);
  r1 = px(i2);
  r2 = py(i2);
  r3 = pz(i2);
  //C
  a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
  b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
  c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
  d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
  ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
  f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
  //C
  //C       make sure particle 2 formed early
  h = a + b;
  if (h > 0e0) {
    g = a;
    a = -b;
    b = -g;
    //C
    g = c;
    c = d;
    d = g;
    //C
    i1 = j;
    i2 = i;
  }
  //C
  //C       check the approaching criteria
  if (allok) {
    //C
    vp = a * d - b * ee;
    //C
    allok = allok && vp < 0e0;
    //C
  }
  //C
  //C       check the closest approach distance criteria
  if (allok) {
    //C
    dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 * a * b *
      ee) / (fem::pow2(ee) - c * d);
    //C
    allok = allok && dm2 < cmn.cutof2;
    //C
  }
  //C
  //C       check the time criteria
  if (allok) {
    //C
    tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
    tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
    if (iordsc == 20) {
      tm = fem::min(tc1, tc2);
    }
    else if (iordsc == 21) {
      tm = 0.5e0 * (tc1 + tc2);
    }
    else {
      tm = fem::max(tc1, tc2);
    }
    //C
    allok = allok && tm > ft(i) && tm > ft(j);
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (h > 0e0) {
    t1 = tc2;
    t2 = tc1;
  }
  else {
    t1 = tc1;
    t2 = tc2;
  }
  //C
}

struct isco3_save
{
  double dgx;
  double dgy;
  double dgz;
  double dm2;
  double dt;
  double dvx;
  double dvy;
  double dvz;
  double dx;
  double dy;
  double dz;
  double e1;
  double e2;
  int i1;
  int i2;
  double px1;
  double px2;
  double py1;
  double py2;
  double pz1;
  double pz2;
  double rts2;
  double v2;
  double vp;
  double vx1;
  double vy1;
  double vz1;

  isco3_save() :
    dgx(fem::double0),
    dgy(fem::double0),
    dgz(fem::double0),
    dm2(fem::double0),
    dt(fem::double0),
    dvx(fem::double0),
    dvy(fem::double0),
    dvz(fem::double0),
    dx(fem::double0),
    dy(fem::double0),
    dz(fem::double0),
    e1(fem::double0),
    e2(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    px1(fem::double0),
    px2(fem::double0),
    py1(fem::double0),
    py2(fem::double0),
    pz1(fem::double0),
    pz2(fem::double0),
    rts2(fem::double0),
    v2(fem::double0),
    vp(fem::double0),
    vx1(fem::double0),
    vy1(fem::double0),
    vz1(fem::double0)
  {}
};

void
isco3(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco3);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& dgx = sve.dgx;
  double& dgy = sve.dgy;
  double& dgz = sve.dgz;
  double& dm2 = sve.dm2;
  double& dt = sve.dt;
  double& dvx = sve.dvx;
  double& dvy = sve.dvy;
  double& dvz = sve.dvz;
  double& dx = sve.dx;
  double& dy = sve.dy;
  double& dz = sve.dz;
  double& e1 = sve.e1;
  double& e2 = sve.e2;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  double& px1 = sve.px1;
  double& px2 = sve.px2;
  double& py1 = sve.py1;
  double& py2 = sve.py2;
  double& pz1 = sve.pz1;
  double& pz2 = sve.pz2;
  double& rts2 = sve.rts2;
  double& v2 = sve.v2;
  double& vp = sve.vp;
  double& vx1 = sve.vx1;
  double& vy1 = sve.vy1;
  double& vz1 = sve.vz1;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  if (ft(i) >= ft(j)) {
    i1 = j;
    i2 = i;
  }
  else {
    i1 = i;
    i2 = j;
  }
  //C
  if (allok) {
    //C
    t1 = ft(i1);
    vx1 = vx(i1);
    vy1 = vy(i1);
    vz1 = vz(i1);
    //C
    t2 = ft(i2);
    //C
    dvx = vx(i2) - vx1;
    dvy = vy(i2) - vy1;
    dvz = vz(i2) - vz1;
    //C
    dt = t2 - t1;
    //C
    dx = gx(i2) - gx(i1) - vx1 * dt;
    dy = gy(i2) - gy(i1) - vy1 * dt;
    dz = gz(i2) - gz(i1) - vz1 * dt;
    //C
    vp = dvx * dx + dvy * dy + dvz * dz;
    //C
    allok = allok && vp < 0e0;
    //C
  }
  //C
  if (allok) {
    //C
    v2 = dvx * dvx + dvy * dvy + dvz * dvz;
    //C
    if (v2 == 0e0) {
      tm = tlarge;
    }
    else {
      tm = t2 - vp / v2;
    }
    //C
    //C       note now tm is the absolute time
    //C
    allok = allok && tm > t1 && tm > t2;
    //C
  }
  //C
  if (allok) {
    //C
    dgx = dx - dvx * t2;
    dgy = dy - dvy * t2;
    dgz = dz - dvz * t2;
    //C
    dm2 = -v2 * fem::pow2(tm) + dgx * dgx + dgy * dgy + dgz * dgz;
    //C
    allok = allok && dm2 < cmn.cutof2;
    //C
  }
  //C
  if (allok) {
    //C
    e1 = e(i1);
    px1 = px(i1);
    py1 = py(i1);
    pz1 = pz(i1);
    e2 = e(i2);
    px2 = px(i2);
    py2 = py(i2);
    pz2 = pz(i2);
    //C
    rts2 = fem::pow2((e1 + e2)) - fem::pow2((px1 + px2)) - fem::pow2((
      py1 + py2)) - fem::pow2((pz1 + pz2));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco4_save
{
  double a;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  int i1;
  int i2;
  int icels1;
  int icels2;
  int ii1;
  int ii2;
  int jj1;
  int jj2;
  int kk1;
  int kk2;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc2;
  double vp;

  isco4_save() :
    a(fem::double0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    icels1(fem::int0),
    icels2(fem::int0),
    ii1(fem::int0),
    ii2(fem::int0),
    jj1(fem::int0),
    jj2(fem::int0),
    kk1(fem::int0),
    kk2(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc2(fem::double0),
    vp(fem::double0)
  {}
};

void
isco4(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco4);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& icels1 = sve.icels1;
  int& icels2 = sve.icels2;
  int& ii1 = sve.ii1;
  int& ii2 = sve.ii2;
  int& jj1 = sve.jj1;
  int& jj2 = sve.jj2;
  int& kk1 = sve.kk1;
  int& kk2 = sve.kk2;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc2 = sve.tc2;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  //C
  icels1 = icels(i);
  ii1 = icels1 / 10000;
  jj1 = (icels1 - ii1 * 10000) / 100;
  kk1 = icels1 - ii1 * 10000 - jj1 * 100;
  icels2 = icels(j);
  ii2 = icels2 / 10000;
  jj2 = (icels2 - ii2 * 10000) / 100;
  kk2 = icels2 - ii2 * 10000 - jj2 * 100;
  //C
  i1 = i;
  i2 = j;
  //C
  p4 = ft(i2) - ft(i1);
  p1 = gx(i2) - gx(i1);
  p2 = gy(i2) - gy(i1);
  p3 = gz(i2) - gz(i1);
  //C
  if (ii1 - ii2 > 5) {
    p1 += 10e0 * size1;
  }
  else if (ii1 - ii2 <  - 5) {
    p1 = p1 - 10e0 * size1;
  }
  if (jj1 - jj2 > 5) {
    p2 += 10e0 * size2;
  }
  else if (jj1 - jj2 <  - 5) {
    p2 = p2 - 10e0 * size2;
  }
  if (kk1 - kk2 > 5) {
    p3 += 10e0 * size3;
  }
  else if (kk1 - kk2 <  - 5) {
    p3 = p3 - 10e0 * size3;
  }
  //C
  q4 = e(i1);
  q1 = px(i1);
  q2 = py(i1);
  q3 = pz(i1);
  //C
  r4 = e(i2);
  r1 = px(i2);
  r2 = py(i2);
  r3 = pz(i2);
  //C
  a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
  b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
  c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
  d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
  ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
  f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
  //C
  //C       make sure particle 2 formed early
  h = a + b;
  if (h > 0e0) {
    g = a;
    a = -b;
    b = -g;
    //C
    g = c;
    c = d;
    d = g;
    //C
    i1 = j;
    i2 = i;
  }
  //C
  //C       check the approaching criteria
  if (allok) {
    //C
    vp = a * d - b * ee;
    //C
    allok = allok && vp < 0e0;
    //C
  }
  //C
  //C       check the closest approach distance criteria
  if (allok) {
    //C
    dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 * a * b *
      ee) / (fem::pow2(ee) - c * d);
    //C
    allok = allok && dm2 < cmn.cutof2;
    //C
  }
  //C
  //C       check the time criteria
  if (allok) {
    //C
    tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
    tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
    tm = 0.5e0 * (tc1 + tc2);
    //C
    allok = allok && tm > ft(i) && tm > ft(j);
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (h > 0e0) {
    t1 = tm;
    t2 = tm;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco5_save
{
  double a;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  int i1;
  int i2;
  int icels1;
  int icels2;
  int ii1;
  int ii2;
  int jj1;
  int jj2;
  int kk1;
  int kk2;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc2;
  double vp;

  isco5_save() :
    a(fem::double0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    icels1(fem::int0),
    icels2(fem::int0),
    ii1(fem::int0),
    ii2(fem::int0),
    jj1(fem::int0),
    jj2(fem::int0),
    kk1(fem::int0),
    kk2(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc2(fem::double0),
    vp(fem::double0)
  {}
};

void
isco5(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco5);
  // COMMON para5
  int& iordsc = cmn.iordsc;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& icels1 = sve.icels1;
  int& icels2 = sve.icels2;
  int& ii1 = sve.ii1;
  int& ii2 = sve.ii2;
  int& jj1 = sve.jj1;
  int& jj2 = sve.jj2;
  int& kk1 = sve.kk1;
  int& kk2 = sve.kk2;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc2 = sve.tc2;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  //C
  icels1 = icels(i);
  ii1 = icels1 / 10000;
  jj1 = (icels1 - ii1 * 10000) / 100;
  kk1 = icels1 - ii1 * 10000 - jj1 * 100;
  icels2 = icels(j);
  ii2 = icels2 / 10000;
  jj2 = (icels2 - ii2 * 10000) / 100;
  kk2 = icels2 - ii2 * 10000 - jj2 * 100;
  //C
  i1 = i;
  i2 = j;
  //C
  p4 = ft(i2) - ft(i1);
  p1 = gx(i2) - gx(i1);
  p2 = gy(i2) - gy(i1);
  p3 = gz(i2) - gz(i1);
  //C
  if (ii1 - ii2 > 5) {
    p1 += 10e0 * size1;
  }
  else if (ii1 - ii2 <  - 5) {
    p1 = p1 - 10e0 * size1;
  }
  if (jj1 - jj2 > 5) {
    p2 += 10e0 * size2;
  }
  else if (jj1 - jj2 <  - 5) {
    p2 = p2 - 10e0 * size2;
  }
  if (kk1 - kk2 > 5) {
    p3 += 10e0 * size3;
  }
  else if (kk1 - kk2 <  - 5) {
    p3 = p3 - 10e0 * size3;
  }
  //C
  q4 = e(i1);
  q1 = px(i1);
  q2 = py(i1);
  q3 = pz(i1);
  //C
  r4 = e(i2);
  r1 = px(i2);
  r2 = py(i2);
  r3 = pz(i2);
  //C
  a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
  b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
  c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
  d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
  ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
  f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
  //C
  //C       make sure particle 2 formed early
  h = a + b;
  if (h > 0e0) {
    g = a;
    a = -b;
    b = -g;
    //C
    g = c;
    c = d;
    d = g;
    //C
    i1 = j;
    i2 = i;
  }
  //C
  //C       check the approaching criteria
  if (allok) {
    //C
    vp = a * d - b * ee;
    //C
    allok = allok && vp < 0e0;
    //C
  }
  //C
  //C       check the closest approach distance criteria
  if (allok) {
    //C
    dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 * a * b *
      ee) / (fem::pow2(ee) - c * d);
    //C
    allok = allok && dm2 < cmn.cutof2;
    //C
  }
  //C
  //C       check the time criteria
  if (allok) {
    //C
    tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
    tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
    if (iordsc == 20) {
      tm = fem::min(tc1, tc2);
    }
    else if (iordsc == 21) {
      tm = 0.5e0 * (tc1 + tc2);
    }
    else {
      tm = fem::max(tc1, tc2);
    }
    //C
    allok = allok && tm > ft(i) && tm > ft(j);
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (h > 0e0) {
    t1 = tc2;
    t2 = tc1;
  }
  else {
    t1 = tc1;
    t2 = tc2;
  }
  //C
}

struct isco6_save
{
  double dgx;
  double dgy;
  double dgz;
  double dm2;
  double dt;
  double dvx;
  double dvy;
  double dvz;
  double dx;
  double dy;
  double dz;
  double e1;
  double e2;
  int i1;
  int i2;
  int icels1;
  int icels2;
  int ii1;
  int ii2;
  int jj1;
  int jj2;
  int kk1;
  int kk2;
  double px1;
  double px2;
  double py1;
  double py2;
  double pz1;
  double pz2;
  double rts2;
  double v2p;
  double vp;
  double vx1;
  double vy1;
  double vz1;

  isco6_save() :
    dgx(fem::double0),
    dgy(fem::double0),
    dgz(fem::double0),
    dm2(fem::double0),
    dt(fem::double0),
    dvx(fem::double0),
    dvy(fem::double0),
    dvz(fem::double0),
    dx(fem::double0),
    dy(fem::double0),
    dz(fem::double0),
    e1(fem::double0),
    e2(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    icels1(fem::int0),
    icels2(fem::int0),
    ii1(fem::int0),
    ii2(fem::int0),
    jj1(fem::int0),
    jj2(fem::int0),
    kk1(fem::int0),
    kk2(fem::int0),
    px1(fem::double0),
    px2(fem::double0),
    py1(fem::double0),
    py2(fem::double0),
    pz1(fem::double0),
    pz2(fem::double0),
    rts2(fem::double0),
    v2p(fem::double0),
    vp(fem::double0),
    vx1(fem::double0),
    vy1(fem::double0),
    vz1(fem::double0)
  {}
};

void
isco6(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco6);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& dgx = sve.dgx;
  double& dgy = sve.dgy;
  double& dgz = sve.dgz;
  double& dm2 = sve.dm2;
  double& dt = sve.dt;
  double& dvx = sve.dvx;
  double& dvy = sve.dvy;
  double& dvz = sve.dvz;
  double& dx = sve.dx;
  double& dy = sve.dy;
  double& dz = sve.dz;
  double& e1 = sve.e1;
  double& e2 = sve.e2;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& icels1 = sve.icels1;
  int& icels2 = sve.icels2;
  int& ii1 = sve.ii1;
  int& ii2 = sve.ii2;
  int& jj1 = sve.jj1;
  int& jj2 = sve.jj2;
  int& kk1 = sve.kk1;
  int& kk2 = sve.kk2;
  double& px1 = sve.px1;
  double& px2 = sve.px2;
  double& py1 = sve.py1;
  double& py2 = sve.py2;
  double& pz1 = sve.pz1;
  double& pz2 = sve.pz2;
  double& rts2 = sve.rts2;
  double& v2p = sve.v2p;
  double& vp = sve.vp;
  double& vx1 = sve.vx1;
  double& vy1 = sve.vy1;
  double& vz1 = sve.vz1;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  if (ft(i) >= ft(j)) {
    i1 = j;
    i2 = i;
  }
  else {
    i1 = i;
    i2 = j;
  }
  //C
  icels1 = icels(i1);
  ii1 = icels1 / 10000;
  jj1 = (icels1 - ii1 * 10000) / 100;
  kk1 = icels1 - ii1 * 10000 - jj1 * 100;
  icels2 = icels(i2);
  ii2 = icels2 / 10000;
  jj2 = (icels2 - ii2 * 10000) / 100;
  kk2 = icels2 - ii2 * 10000 - jj2 * 100;
  //C
  if (allok) {
    //C
    t1 = ft(i1);
    vx1 = vx(i1);
    vy1 = vy(i1);
    vz1 = vz(i1);
    //C
    t2 = ft(i2);
    //C
    dvx = vx(i2) - vx1;
    dvy = vy(i2) - vy1;
    dvz = vz(i2) - vz1;
    //C
    dt = t2 - t1;
    //C
    dx = gx(i2) - gx(i1) - vx1 * dt;
    dy = gy(i2) - gy(i1) - vy1 * dt;
    dz = gz(i2) - gz(i1) - vz1 * dt;
    //C
    if (ii1 - ii2 > 5) {
      dx += 10e0 * size1;
    }
    else if (ii1 - ii2 <  - 5) {
      dx = dx - 10e0 * size1;
    }
    //C
    if (jj1 - jj2 > 5) {
      dy += 10e0 * size2;
    }
    else if (jj1 - jj2 <  - 5) {
      dy = dy - 10e0 * size2;
    }
    //C
    if (kk1 - kk2 > 5) {
      dz += 10e0 * size3;
    }
    else if (kk1 - kk2 <  - 5) {
      dz = dz - 10e0 * size3;
    }
    //C
    vp = dvx * dx + dvy * dy + dvz * dz;
    //C
    allok = allok && vp < 0e0;
    //C
  }
  //C
  if (allok) {
    //C
    v2p = dvx * dvx + dvy * dvy + dvz * dvz;
    //C
    if (v2p == 0e0) {
      tm = tlarge;
    }
    else {
      tm = t2 - vp / v2p;
    }
    //C
    //C       note now tm is the absolute time
    //C
    allok = allok && tm > t1 && tm > t2;
    //C
  }
  //C
  if (allok) {
    //C
    dgx = dx - dvx * t2;
    dgy = dy - dvy * t2;
    dgz = dz - dvz * t2;
    //C
    dm2 = -v2p * fem::pow2(tm) + dgx * dgx + dgy * dgy + dgz * dgz;
    //C
    allok = allok && dm2 < cmn.cutof2;
    //C
  }
  //C
  if (allok) {
    //C
    e1 = e(i1);
    px1 = px(i1);
    py1 = py(i1);
    pz1 = pz(i1);
    e2 = e(i2);
    px2 = px(i2);
    py2 = py(i2);
    pz2 = pz(i2);
    //C
    rts2 = fem::pow2((e1 + e2)) - fem::pow2((px1 + px2)) - fem::pow2((
      py1 + py2)) - fem::pow2((pz1 + pz2));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco7_save
{
  double a;
  bool allokp;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  int i1;
  int i2;
  int ii;
  int jj;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc2;
  double tmp;
  double vp;

  isco7_save() :
    a(fem::double0),
    allokp(fem::bool0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    ii(fem::int0),
    jj(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc2(fem::double0),
    tmp(fem::double0),
    vp(fem::double0)
  {}
};

void
isco7(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco7);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  bool& allokp = sve.allokp;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& ii = sve.ii;
  int& jj = sve.jj;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc2 = sve.tc2;
  double& tmp = sve.tmp;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  //C
  tm = tlarge;
  //C
  if (allok) {
    FEM_DO_SAFE(ii, -1, 1) {
      FEM_DO_SAFE(jj, -1, 1) {
        //C
        allokp = true;
        //C
        i1 = i;
        i2 = j;
        //C
        p4 = ft(j) - ft(i);
        p1 = gx(j) - gx(i);
        p2 = gy(j) - gy(i);
        p3 = gz(j) - gz(i);
        //C
        p1 += ii * 10e0 * cmn.size1;
        p2 += jj * 10e0 * cmn.size2;
        //C
        q4 = e(i);
        q1 = px(i);
        q2 = py(i);
        q3 = pz(i);
        //C
        r4 = e(j);
        r1 = px(j);
        r2 = py(j);
        r3 = pz(j);
        //C
        a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
        b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
        c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
        d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
        ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
        f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
        //C
        //C       make sure particle 2 formed early
        h = a + b;
        if (h > 0e0) {
          g = a;
          a = -b;
          b = -g;
          g = c;
          c = d;
          d = g;
          i1 = j;
          i2 = i;
        }
        //C
        //C       check the approaching criteria
        if (allokp) {
          vp = a * d - b * ee;
          allokp = allokp && vp < 0e0;
        }
        //C
        //C       check the closest approach distance criteria
        if (allokp) {
          dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 * a *
            b * ee) / (fem::pow2(ee) - c * d);
          allokp = allokp && dm2 < cmn.cutof2;
        }
        //C
        //C       check the time criteria
        if (allokp) {
          tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
          tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
          tmp = 0.5e0 * (tc1 + tc2);
          allokp = allokp && tmp > ft(i) && tmp > ft(j);
        }
        //C
        if (allokp && tmp < tm) {
          tm = tmp;
          cmn.jxa = ii;
          cmn.jya = jj;
          //Cd                    dgxa(j) = ii * 10d0 * size1
          //Cd                    dgya(j) = jj * 10d0 * size2
          //Cd                    dgxa(i) = - dgxa(j)
          //Cd                    dgya(i) = - dgya(j)
        }
        //C
      }
    }
    //C
    if (tm == tlarge) {
      allok = false;
    }
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    q4 = e(i1);
    q1 = px(i1);
    q2 = py(i1);
    q3 = pz(i1);
    //C
    r4 = e(i2);
    r1 = px(i2);
    r2 = py(i2);
    r3 = pz(i2);
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (h > 0e0) {
    t1 = tm;
    t2 = tm;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco8_save
{
  double a;
  bool allokp;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  double ha;
  int i1;
  int i2;
  int ii;
  int jj;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc1a;
  double tc2;
  double tc2a;
  double tmp;
  double vp;

  isco8_save() :
    a(fem::double0),
    allokp(fem::bool0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    ha(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    ii(fem::int0),
    jj(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc1a(fem::double0),
    tc2(fem::double0),
    tc2a(fem::double0),
    tmp(fem::double0),
    vp(fem::double0)
  {}
};

void
isco8(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco8);
  // COMMON para5
  int& iordsc = cmn.iordsc;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  bool& allokp = sve.allokp;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  double& ha = sve.ha;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& ii = sve.ii;
  int& jj = sve.jj;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc1a = sve.tc1a;
  double& tc2 = sve.tc2;
  double& tc2a = sve.tc2a;
  double& tmp = sve.tmp;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  //C
  tm = tlarge;
  //C
  if (allok) {
    FEM_DO_SAFE(ii, -1, 1) {
      FEM_DO_SAFE(jj, -1, 1) {
        //C
        allokp = true;
        //C
        i1 = i;
        i2 = j;
        //C
        p4 = ft(j) - ft(i);
        p1 = gx(j) - gx(i);
        p2 = gy(j) - gy(i);
        p3 = gz(j) - gz(i);
        //C
        p1 += ii * 10e0 * cmn.size1;
        p2 += jj * 10e0 * cmn.size2;
        //C
        q4 = e(i);
        q1 = px(i);
        q2 = py(i);
        q3 = pz(i);
        //C
        r4 = e(j);
        r1 = px(j);
        r2 = py(j);
        r3 = pz(j);
        //C
        a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
        b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
        c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
        d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
        ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
        f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
        //C
        //C       make sure particle 2 formed early
        h = a + b;
        if (h > 0e0) {
          g = a;
          a = -b;
          b = -g;
          g = c;
          c = d;
          d = g;
          i1 = j;
          i2 = i;
        }
        //C
        //C       check the approaching criteria
        if (allokp) {
          vp = a * d - b * ee;
          allokp = allokp && vp < 0e0;
        }
        //C
        //C       check the closest approach distance criteria
        if (allokp) {
          dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 * a *
            b * ee) / (fem::pow2(ee) - c * d);
          allokp = allokp && dm2 < cmn.cutof2;
        }
        //C
        //C       check the time criteria
        if (allokp) {
          tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
          tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
          if (iordsc == 20) {
            tmp = fem::min(tc1, tc2);
          }
          else if (iordsc == 21) {
            tmp = 0.5e0 * (tc1 + tc2);
          }
          else {
            tmp = fem::max(tc1, tc2);
          }
          allokp = allokp && tmp > ft(i) && tmp > ft(j);
        }
        //C
        if (allokp && tmp < tm) {
          tm = tmp;
          cmn.jxa = ii;
          cmn.jya = jj;
          ha = h;
          tc1a = tc1;
          tc2a = tc2;
          //Cd                    dgxa(j) = ii * 10d0 * size1
          //Cd                    dgya(j) = jj * 10d0 * size2
          //Cd                    dgxa(i) = - dgxa(j)
          //Cd                    dgya(i) = - dgya(j)
        }
        //C
      }
    }
    //C
    if (tm == tlarge) {
      allok = false;
    }
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    q4 = e(i1);
    q1 = px(i1);
    q2 = py(i1);
    q3 = pz(i1);
    //C
    r4 = e(i2);
    r1 = px(i2);
    r2 = py(i2);
    r3 = pz(i2);
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (ha > 0e0) {
    t1 = tc2a;
    t2 = tc1a;
  }
  else {
    t1 = tc1a;
    t2 = tc2a;
  }
  //C
}

struct isco9_save
{
  bool allokp;
  double dgx;
  double dgy;
  double dgz;
  double dm2;
  double dt;
  double dvx;
  double dvy;
  double dvz;
  double dx;
  double dy;
  double dz;
  double e1;
  double e2;
  int i1;
  int i2;
  int ii;
  int isign;
  int jj;
  double px1;
  double px2;
  double py1;
  double py2;
  double pz1;
  double pz2;
  double rts2;
  double tmp;
  double vp;
  double vx1;
  double vy1;
  double vz1;

  isco9_save() :
    allokp(fem::bool0),
    dgx(fem::double0),
    dgy(fem::double0),
    dgz(fem::double0),
    dm2(fem::double0),
    dt(fem::double0),
    dvx(fem::double0),
    dvy(fem::double0),
    dvz(fem::double0),
    dx(fem::double0),
    dy(fem::double0),
    dz(fem::double0),
    e1(fem::double0),
    e2(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    ii(fem::int0),
    isign(fem::int0),
    jj(fem::int0),
    px1(fem::double0),
    px2(fem::double0),
    py1(fem::double0),
    py2(fem::double0),
    pz1(fem::double0),
    pz2(fem::double0),
    rts2(fem::double0),
    tmp(fem::double0),
    vp(fem::double0),
    vx1(fem::double0),
    vy1(fem::double0),
    vz1(fem::double0)
  {}
};

void
isco9(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco9);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist3
  double& v2 = cmn.v2;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  bool& allokp = sve.allokp;
  double& dgx = sve.dgx;
  double& dgy = sve.dgy;
  double& dgz = sve.dgz;
  double& dm2 = sve.dm2;
  double& dt = sve.dt;
  double& dvx = sve.dvx;
  double& dvy = sve.dvy;
  double& dvz = sve.dvz;
  double& dx = sve.dx;
  double& dy = sve.dy;
  double& dz = sve.dz;
  double& e1 = sve.e1;
  double& e2 = sve.e2;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& ii = sve.ii;
  int& isign = sve.isign;
  int& jj = sve.jj;
  double& px1 = sve.px1;
  double& px2 = sve.px2;
  double& py1 = sve.py1;
  double& py2 = sve.py2;
  double& pz1 = sve.pz1;
  double& pz2 = sve.pz2;
  double& rts2 = sve.rts2;
  double& tmp = sve.tmp;
  double& vp = sve.vp;
  double& vx1 = sve.vx1;
  double& vy1 = sve.vy1;
  double& vz1 = sve.vz1;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  if (ft(i) >= ft(j)) {
    i1 = j;
    i2 = i;
    isign = -1;
  }
  else {
    i1 = i;
    i2 = j;
    isign = 1;
  }
  //C
  if (allok) {
    tm = tlarge;
    //C
    t1 = ft(i1);
    vx1 = vx(i1);
    vy1 = vy(i1);
    vz1 = vz(i1);
    //C
    t2 = ft(i2);
    //C
    dvx = vx(i2) - vx1;
    dvy = vy(i2) - vy1;
    dvz = vz(i2) - vz1;
    //C
    dt = t2 - t1;
    //C
    FEM_DO_SAFE(ii, -1, 1) {
      FEM_DO_SAFE(jj, -1, 1) {
        //C
        allokp = true;
        //C
        dx = gx(i2) - gx(i1) - vx1 * dt;
        dy = gy(i2) - gy(i1) - vy1 * dt;
        dz = gz(i2) - gz(i1) - vz1 * dt;
        //C
        dx += ii * 10e0 * cmn.size1;
        dy += jj * 10e0 * cmn.size2;
        //C
        vp = dvx * dx + dvy * dy + dvz * dz;
        //C
        allokp = allokp && vp < 0e0;
        //C
        if (allokp) {
          //C
          v2 = dvx * dvx + dvy * dvy + dvz * dvz;
          //C
          if (v2 == 0e0) {
            tmp = tlarge;
          }
          else {
            tmp = t2 - vp / v2;
          }
          //C
          //C       note now tm is the absolute time
          //C
          allokp = allokp && tmp > t1 && tmp > t2;
          //C
        }
        //C
        if (allokp) {
          //C
          dgx = dx - dvx * t2;
          dgy = dy - dvy * t2;
          dgz = dz - dvz * t2;
          //C
          dm2 = -v2 * fem::pow2(tmp) + dgx * dgx + dgy * dgy + dgz * dgz;
          //C
          allokp = allokp && dm2 < cmn.cutof2;
          //C
        }
        //C
        if (allokp && tmp < tm) {
          tm = tmp;
          cmn.jxa = isign * ii;
          cmn.jya = isign * jj;
        }
        //C
      }
    }
    //C
    if (tm == tlarge) {
      allok = false;
    }
  }
  //C
  if (allok) {
    //C
    e1 = e(i1);
    px1 = px(i1);
    py1 = py(i1);
    pz1 = pz(i1);
    e2 = e(i2);
    px2 = px(i2);
    py2 = py(i2);
    pz2 = pz(i2);
    //C
    rts2 = fem::pow2((e1 + e2)) - fem::pow2((px1 + px2)) - fem::pow2((
      py1 + py2)) - fem::pow2((pz1 + pz2));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco10_save
{
  double a;
  bool allokp;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  int i1;
  int i2;
  int ii;
  int jj;
  int kk;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc2;
  double tmp;
  double vp;

  isco10_save() :
    a(fem::double0),
    allokp(fem::bool0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    ii(fem::int0),
    jj(fem::int0),
    kk(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc2(fem::double0),
    tmp(fem::double0),
    vp(fem::double0)
  {}
};

void
isco10(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco10);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  bool& allokp = sve.allokp;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& ii = sve.ii;
  int& jj = sve.jj;
  int& kk = sve.kk;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc2 = sve.tc2;
  double& tmp = sve.tmp;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  //C
  tm = tlarge;
  //C
  if (allok) {
    FEM_DO_SAFE(ii, -1, 1) {
      FEM_DO_SAFE(jj, -1, 1) {
        FEM_DO_SAFE(kk, -1, 1) {
          allokp = true;
          //C
          i1 = i;
          i2 = j;
          //C
          p4 = ft(j) - ft(i);
          p1 = gx(j) - gx(i);
          p2 = gy(j) - gy(i);
          p3 = gz(j) - gz(i);
          //C
          p1 += ii * 10e0 * cmn.size1;
          p2 += jj * 10e0 * cmn.size2;
          p3 += kk * 10e0 * cmn.size3;
          //C
          q4 = e(i);
          q1 = px(i);
          q2 = py(i);
          q3 = pz(i);
          //C
          r4 = e(j);
          r1 = px(j);
          r2 = py(j);
          r3 = pz(j);
          //C
          a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
          b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
          c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
          d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
          ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
          f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
          //C
          //C       make sure particle 2 formed early
          h = a + b;
          if (h > 0e0) {
            g = a;
            a = -b;
            b = -g;
            g = c;
            c = d;
            d = g;
            i1 = j;
            i2 = i;
          }
          //C
          //C       check the approaching criteria
          if (allokp) {
            vp = a * d - b * ee;
            allokp = allokp && vp < 0e0;
          }
          //C
          //C       check the closest approach distance criteria
          if (allokp) {
            dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 *
              a * b * ee) / (fem::pow2(ee) - c * d);
            allokp = allokp && dm2 < cmn.cutof2;
          }
          //C
          //C       check the time criteria
          if (allokp) {
            tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
            tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
            tmp = 0.5e0 * (tc1 + tc2);
            allokp = allokp && tmp > ft(i) && tmp > ft(j);
          }
          //C
          if (allokp && tmp < tm) {
            tm = tmp;
            cmn.jxa = ii;
            cmn.jya = jj;
            cmn.jza = kk;
            //Cd                    dgxa(j) = ii * 10d0 * size1
            //Cd                    dgya(j) = jj * 10d0 * size2
            //Cd                    dgxa(i) = - dgxa(j)
            //Cd                    dgya(i) = - dgya(j)
          }
          //C
        }
      }
    }
    //C
    if (tm == tlarge) {
      allok = false;
    }
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    q4 = e(i1);
    q1 = px(i1);
    q2 = py(i1);
    q3 = pz(i1);
    //C
    r4 = e(i2);
    r1 = px(i2);
    r2 = py(i2);
    r3 = pz(i2);
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (h > 0e0) {
    t1 = tm;
    t2 = tm;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco11_save
{
  double a;
  bool allokp;
  double b;
  double c;
  double d;
  double dm2;
  double ee;
  double f;
  double g;
  double h;
  double ha;
  int i1;
  int i2;
  int ii;
  int jj;
  int kk;
  double p1;
  double p2;
  double p3;
  double p4;
  double q1;
  double q2;
  double q3;
  double q4;
  double r1;
  double r2;
  double r3;
  double r4;
  double rts2;
  double tc1;
  double tc1a;
  double tc2;
  double tc2a;
  double tmp;
  double vp;

  isco11_save() :
    a(fem::double0),
    allokp(fem::bool0),
    b(fem::double0),
    c(fem::double0),
    d(fem::double0),
    dm2(fem::double0),
    ee(fem::double0),
    f(fem::double0),
    g(fem::double0),
    h(fem::double0),
    ha(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    ii(fem::int0),
    jj(fem::int0),
    kk(fem::int0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    p4(fem::double0),
    q1(fem::double0),
    q2(fem::double0),
    q3(fem::double0),
    q4(fem::double0),
    r1(fem::double0),
    r2(fem::double0),
    r3(fem::double0),
    r4(fem::double0),
    rts2(fem::double0),
    tc1(fem::double0),
    tc1a(fem::double0),
    tc2(fem::double0),
    tc2a(fem::double0),
    tmp(fem::double0),
    vp(fem::double0)
  {}
};

void
isco11(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco11);
  // COMMON para5
  int& iordsc = cmn.iordsc;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  double& a = sve.a;
  bool& allokp = sve.allokp;
  double& b = sve.b;
  double& c = sve.c;
  double& d = sve.d;
  double& dm2 = sve.dm2;
  double& ee = sve.ee;
  double& f = sve.f;
  double& g = sve.g;
  double& h = sve.h;
  double& ha = sve.ha;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& ii = sve.ii;
  int& jj = sve.jj;
  int& kk = sve.kk;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& p4 = sve.p4;
  double& q1 = sve.q1;
  double& q2 = sve.q2;
  double& q3 = sve.q3;
  double& q4 = sve.q4;
  double& r1 = sve.r1;
  double& r2 = sve.r2;
  double& r3 = sve.r3;
  double& r4 = sve.r4;
  double& rts2 = sve.rts2;
  double& tc1 = sve.tc1;
  double& tc1a = sve.tc1a;
  double& tc2 = sve.tc2;
  double& tc2a = sve.tc2a;
  double& tmp = sve.tmp;
  double& vp = sve.vp;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  //C       set up numbers for later calculations
  //C
  tm = tlarge;
  //C
  if (allok) {
    FEM_DO_SAFE(ii, -1, 1) {
      FEM_DO_SAFE(jj, -1, 1) {
        FEM_DO_SAFE(kk, -1, 1) {
          //C
          allokp = true;
          //C
          i1 = i;
          i2 = j;
          //C
          p4 = ft(j) - ft(i);
          p1 = gx(j) - gx(i);
          p2 = gy(j) - gy(i);
          p3 = gz(j) - gz(i);
          //C
          p1 += ii * 10e0 * cmn.size1;
          p2 += jj * 10e0 * cmn.size2;
          p3 += kk * 10e0 * cmn.size3;
          //C
          q4 = e(i);
          q1 = px(i);
          q2 = py(i);
          q3 = pz(i);
          //C
          r4 = e(j);
          r1 = px(j);
          r2 = py(j);
          r3 = pz(j);
          //C
          a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
          b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
          c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
          d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
          ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
          f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;
          //C
          //C       make sure particle 2 formed early
          h = a + b;
          if (h > 0e0) {
            g = a;
            a = -b;
            b = -g;
            g = c;
            c = d;
            d = g;
            i1 = j;
            i2 = i;
          }
          //C
          //C       check the approaching criteria
          if (allokp) {
            vp = a * d - b * ee;
            allokp = allokp && vp < 0e0;
          }
          //C
          //C       check the closest approach distance criteria
          if (allokp) {
            dm2 = -f - (fem::pow2(a) * d + fem::pow2(b) * c - 2e0 *
              a * b * ee) / (fem::pow2(ee) - c * d);
            allokp = allokp && dm2 < cmn.cutof2;
          }
          //C
          //C       check the time criteria
          if (allokp) {
            tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (fem::pow2(ee) - c * d);
            tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (fem::pow2(ee) - c * d);
            if (iordsc == 20) {
              tmp = fem::min(tc1, tc2);
            }
            else if (iordsc == 21) {
              tmp = 0.5e0 * (tc1 + tc2);
            }
            else {
              tmp = fem::max(tc1, tc2);
            }
            allokp = allokp && tmp > ft(i) && tmp > ft(j);
          }
          //C
          if (allokp && tmp < tm) {
            tm = tmp;
            cmn.jxa = ii;
            cmn.jya = jj;
            cmn.jza = kk;
            ha = h;
            tc1a = tc1;
            tc2a = tc2;
            //Cd                    dgxa(j) = ii * 10d0 * size1
            //Cd                    dgya(j) = jj * 10d0 * size2
            //Cd                    dgxa(i) = - dgxa(j)
            //Cd                    dgya(i) = - dgya(j)
          }
          //C
        }
      }
    }
    //C
    if (tm == tlarge) {
      allok = false;
    }
    //C
  }
  //C
  //C        check rts cut
  if (allok) {
    //C
    q4 = e(i1);
    q1 = px(i1);
    q2 = py(i1);
    q3 = pz(i1);
    //C
    r4 = e(i2);
    r1 = px(i2);
    r2 = py(i2);
    r3 = pz(i2);
    //C
    rts2 = fem::pow2((q4 + r4)) - fem::pow2((q1 + r1)) - fem::pow2((
      q2 + r2)) - fem::pow2((q3 + r3));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else if (ha > 0e0) {
    t1 = tc2a;
    t2 = tc1a;
  }
  else {
    t1 = tc1a;
    t2 = tc2a;
  }
  //C
}

struct isco12_save
{
  bool allokp;
  double dgx;
  double dgy;
  double dgz;
  double dm2;
  double dt;
  double dvx;
  double dvy;
  double dvz;
  double dx;
  double dy;
  double dz;
  double e1;
  double e2;
  int i1;
  int i2;
  int ii;
  int isign;
  int jj;
  int kk;
  double px1;
  double px2;
  double py1;
  double py2;
  double pz1;
  double pz2;
  double rts2;
  double tmp;
  double vp;
  double vx1;
  double vy1;
  double vz1;

  isco12_save() :
    allokp(fem::bool0),
    dgx(fem::double0),
    dgy(fem::double0),
    dgz(fem::double0),
    dm2(fem::double0),
    dt(fem::double0),
    dvx(fem::double0),
    dvy(fem::double0),
    dvz(fem::double0),
    dx(fem::double0),
    dy(fem::double0),
    dz(fem::double0),
    e1(fem::double0),
    e2(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    ii(fem::int0),
    isign(fem::int0),
    jj(fem::int0),
    kk(fem::int0),
    px1(fem::double0),
    px2(fem::double0),
    py1(fem::double0),
    py2(fem::double0),
    pz1(fem::double0),
    pz2(fem::double0),
    rts2(fem::double0),
    tmp(fem::double0),
    vp(fem::double0),
    vx1(fem::double0),
    vy1(fem::double0),
    vz1(fem::double0)
  {}
};

void
isco12(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco12);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> last(cmn.last, dimension(maxptn));
  // COMMON ilist3
  double& v2 = cmn.v2;
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  bool& allokp = sve.allokp;
  double& dgx = sve.dgx;
  double& dgy = sve.dgy;
  double& dgz = sve.dgz;
  double& dm2 = sve.dm2;
  double& dt = sve.dt;
  double& dvx = sve.dvx;
  double& dvy = sve.dvy;
  double& dvz = sve.dvz;
  double& dx = sve.dx;
  double& dy = sve.dy;
  double& dz = sve.dz;
  double& e1 = sve.e1;
  double& e2 = sve.e2;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& ii = sve.ii;
  int& isign = sve.isign;
  int& jj = sve.jj;
  int& kk = sve.kk;
  double& px1 = sve.px1;
  double& px2 = sve.px2;
  double& py1 = sve.py1;
  double& py2 = sve.py2;
  double& pz1 = sve.pz1;
  double& pz2 = sve.pz2;
  double& rts2 = sve.rts2;
  double& tmp = sve.tmp;
  double& vp = sve.vp;
  double& vx1 = sve.vx1;
  double& vy1 = sve.vy1;
  double& vz1 = sve.vz1;
  //
  //C       this subroutine is used to decide whether there is a collision between
  //C       particle i and j, if there is one allok=1, and tm gives the
  //C       collision time, t1 the collision time for i,
  //C       t2 the collision time for j
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  //C       preventing consecutive collisions
  allok = last(i) != j || last(j) != i;
  //C
  if (ft(i) >= ft(j)) {
    i1 = j;
    i2 = i;
    isign = -1;
  }
  else {
    i1 = i;
    i2 = j;
    isign = 1;
  }
  //C
  if (allok) {
    tm = tlarge;
    //C
    t1 = ft(i1);
    vx1 = vx(i1);
    vy1 = vy(i1);
    vz1 = vz(i1);
    //C
    t2 = ft(i2);
    //C
    dvx = vx(i2) - vx1;
    dvy = vy(i2) - vy1;
    dvz = vz(i2) - vz1;
    //C
    dt = t2 - t1;
    //C
    FEM_DO_SAFE(ii, -1, 1) {
      FEM_DO_SAFE(jj, -1, 1) {
        FEM_DO_SAFE(kk, -1, 1) {
          //C
          allokp = true;
          //C
          dx = gx(i2) - gx(i1) - vx1 * dt;
          dy = gy(i2) - gy(i1) - vy1 * dt;
          dz = gz(i2) - gz(i1) - vz1 * dt;
          //C
          dx += ii * 10e0 * cmn.size1;
          dy += jj * 10e0 * cmn.size2;
          dz += kk * 10e0 * cmn.size3;
          //C
          vp = dvx * dx + dvy * dy + dvz * dz;
          //C
          allokp = allokp && vp < 0e0;
          //C
          if (allokp) {
            //C
            v2 = dvx * dvx + dvy * dvy + dvz * dvz;
            //C
            if (v2 == 0e0) {
              tmp = tlarge;
            }
            else {
              tmp = t2 - vp / v2;
            }
            //C
            //C       note now tm is the absolute time
            //C
            allokp = allokp && tmp > t1 && tmp > t2;
            //C
          }
          //C
          if (allokp) {
            //C
            dgx = dx - dvx * t2;
            dgy = dy - dvy * t2;
            dgz = dz - dvz * t2;
            //C
            dm2 = -v2 * fem::pow2(tmp) + dgx * dgx + dgy * dgy + dgz * dgz;
            //C
            allokp = allokp && dm2 < cmn.cutof2;
            //C
          }
          //C
          if (allokp && tmp < tm) {
            tm = tmp;
            cmn.jxa = isign * ii;
            cmn.jya = isign * jj;
            cmn.jza = isign * kk;
          }
          //C
        }
      }
    }
    //C
    if (tm == tlarge) {
      allok = false;
    }
  }
  //C
  if (allok) {
    //C
    e1 = e(i1);
    px1 = px(i1);
    py1 = py(i1);
    pz1 = pz(i1);
    e2 = e(i2);
    px2 = px(i2);
    py2 = py(i2);
    pz2 = pz(i2);
    //C
    rts2 = fem::pow2((e1 + e2)) - fem::pow2((px1 + px2)) - fem::pow2((
      py1 + py2)) - fem::pow2((pz1 + pz2));
    //C
    allok = allok && rts2 > cmn.rscut2;
  }
  //C
  if (!allok) {
    tm = tlarge;
    t1 = tlarge;
    t2 = tlarge;
  }
  else {
    t1 = tm;
    t2 = tm;
  }
  //C
}

struct isco_save
{
  int iorder;

  isco_save() :
    iorder(fem::int0)
  {}
};

void
isco(
  common& cmn,
  int const& i,
  int const& j,
  bool& allok,
  double& tm,
  double& t1,
  double& t2)
{
  FEM_CMN_SVE(isco);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  //
  // SAVE
  int& iorder = sve.iorder;
  //
  //C
  //Cc      SAVE /para5/
  //C
  iorder = cmn.iordsc / 10;
  if (iconfg == 1) {
    if (iorder == 1) {
      isco1(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 2) {
      isco2(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 3) {
      isco3(cmn, i, j, allok, tm, t1, t2);
    }
  }
  else if (iconfg == 2 || iconfg == 4) {
    if (iorder == 1) {
      isco4(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 2) {
      isco5(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 3) {
      isco6(cmn, i, j, allok, tm, t1, t2);
    }
  }
  else if (iconfg == 3) {
    if (iorder == 1) {
      isco7(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 2) {
      isco8(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 3) {
      isco9(cmn, i, j, allok, tm, t1, t2);
    }
  }
  else if (iconfg == 5) {
    if (iorder == 1) {
      isco10(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 2) {
      isco11(cmn, i, j, allok, tm, t1, t2);
    }
    else if (iorder == 3) {
      isco12(cmn, i, j, allok, tm, t1, t2);
    }
  }
  //C
}

struct mintm_save
{
  bool allok;
  double t1;
  double t2;
  double tm;

  mintm_save() :
    allok(fem::bool0),
    t1(fem::double0),
    t2(fem::double0),
    tm(fem::double0)
  {}
};

void
mintm(
  common& cmn,
  int const& i,
  int const& j,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(mintm);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON aurec2
  const int maxptn = 400001;
  arr_ref<double> dgxa(cmn.dgxa, dimension(maxptn));
  arr_ref<double> dgya(cmn.dgya, dimension(maxptn));
  arr_ref<double> dgza(cmn.dgza, dimension(maxptn));
  // COMMON ilist5
  arr_ref<double> ct(cmn.ct, dimension(maxptn));
  //
  // SAVE
  bool& allok = sve.allok;
  double& t1 = sve.t1;
  double& tm = sve.tm;
  //
  //C       this subroutine is used to check whether particle j has smaller
  //C       collision time with particle i than other particles
  //C       or in other words, update next(i)
  //C
  //C       input i,j,tmin,nc
  //C       output tmin,nc
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  isco(cmn, i, j, allok, tm, t1, sve.t2);
  //C
  if (allok && tm < tmin) {
    tmin = tm;
    ct(i) = t1;
    nc = j;
    if (iconfg == 3 || iconfg == 5) {
      dgxa(i) = -cmn.jxa * 10e0 * cmn.size1;
      dgya(i) = -cmn.jya * 10e0 * cmn.size2;
      if (iconfg == 5) {
        dgza(i) = -cmn.jza * 10e0 * cmn.size3;
      }
    }
  }
  //C
}

struct chcell_save
{
  int j;
  int jj;
  int jmintm;
  int l;

  chcell_save() :
    j(fem::int0),
    jj(fem::int0),
    jmintm(fem::int0),
    l(fem::int0)
  {}
};

void
chcell(
  common& cmn,
  int const& il,
  int const& i1,
  int const& i2,
  int const& i3,
  int const& last0,
  double const& /* t */,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chcell);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON ilist1
  const int maxptn = 400001;
  arr_cref<int> nic(cmn.nic, dimension(maxptn));
  // COMMON ilist2
  arr_cref<int, 3> icel(cmn.icel, dimension(10, 10, 10));
  //
  // SAVE
  int& j = sve.j;
  int& jj = sve.jj;
  int& jmintm = sve.jmintm;
  int& l = sve.l;
  //
  //C       this program is used to check through all the particles, except last0
  //C       in the cell (i1,i2,i3) and see if we can get a particle collision
  //C       with time less than the original input tmin ( the collision time of
  //C       il with the wall
  //C       last0 cas be set to 0 if we don't want to exclude last0
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist4/
  //C
  if (iconfg == 3 || iconfg == 5) {
    jj = cmn.ichkpt;
    FEM_DO_SAFE(j, 1, jj) {
      //C     10/24/02 get rid of argument usage mismatch in mintm():
      jmintm = j;
      if (j != il && j != last0) {
        mintm(cmn, il, jmintm, tmin, nc);
      }
      //C     &          call mintm(il, j, tmin, nc)
      //C
    }
    return;
  }
  //C
  //C       set l
  if (i1 == 11 && i2 == 11 && i3 == 11) {
    l = cmn.icell;
  }
  else {
    l = icel(i1, i2, i3);
  }
  //C
  if (l == 0) {
    return;
  }
  //C
  j = nic(l);
  //C
  //C       if there is only one particle
  if (j == 0) {
    //C
    //C       if it's not il or last0,when last is not wall
    if (l == il || l == last0) {
      return;
    }
    mintm(cmn, il, l, tmin, nc);
    //C
    //C       if there are many particles
  }
  else {
    if (l != il && l != last0) {
      mintm(cmn, il, l, tmin, nc);
    }
    while (j != l) {
      if (j != il && j != last0) {
        mintm(cmn, il, j, tmin, nc);
      }
      j = nic(j);
    }
  }
  //C
}

struct chout_save
{
  int i;
  int j;
  int k;
  int m1;
  int m2;
  int m3;

  chout_save() :
    i(fem::int0),
    j(fem::int0),
    k(fem::int0),
    m1(fem::int0),
    m2(fem::int0),
    m3(fem::int0)
  {}
};

void
chout(
  common& cmn,
  int const& l,
  int const& last0,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chout);
  // SAVE
  int& i = sve.i;
  int& j = sve.j;
  int& k = sve.k;
  int& m1 = sve.m1;
  int& m2 = sve.m2;
  int& m3 = sve.m3;
  //
  //C       this subroutine is used to check the surface when the colliding
  //C       particle is outside the cube
  //C
  //Cc      SAVE /prec2/
  //C
  m1 = 11;
  m2 = 11;
  m3 = 11;
  chcell(cmn, l, m1, m2, m3, last0, t, tmin, nc);
  //C
  FEM_DO_SAFE(i, 1, 10) {
    FEM_DO_SAFE(j, 1, 10) {
      FEM_DO_SAFE(k, 1, 10) {
        if (i == 1 || i == 10 || j == 1 || j == 10 || k == 1 || k == 10) {
          chcell(cmn, l, i, j, k, last0, t, tmin, nc);
        }
      }
    }
  }
  //C
}

struct chin1_save
{
  int i;
  int itest;
  int j;
  int k;
  int m1;
  int m2;
  int m3;

  chin1_save() :
    i(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0),
    m1(fem::int0),
    m2(fem::int0),
    m3(fem::int0)
  {}
};

void
chin1(
  common& cmn,
  int const& l,
  int const& i1,
  int const& i2,
  int const& i3,
  int const& last0,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chin1);
  // SAVE
  int& i = sve.i;
  int& itest = sve.itest;
  int& j = sve.j;
  int& k = sve.k;
  int& m1 = sve.m1;
  int& m2 = sve.m2;
  int& m3 = sve.m3;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube
  //C       and update the afftected particles through chcell
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0;
  //C
  FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i >= 1 && i <= 10 && j >= 1 && j <= 10 && k >= 1 && k <= 10) {
          chcell(cmn, l, i, j, k, last0, t, tmin, nc);
        }
        else if (itest == 0) {
          m1 = 11;
          m2 = 11;
          m3 = 11;
          chcell(cmn, l, m1, m2, m3, last0, t, tmin, nc);
          itest = 1;
        }
      }
    }
  }
  //C
}

struct chin2_save
{
  int i;
  int ia;
  int ib;
  int ic;
  int itest;
  int j;
  int k;

  chin2_save() :
    i(fem::int0),
    ia(fem::int0),
    ib(fem::int0),
    ic(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
chin2(
  common& cmn,
  int const& l,
  int const& i1,
  int const& i2,
  int const& i3,
  int const& last0,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chin2);
  // SAVE
  int& i = sve.i;
  int& ia = sve.ia;
  int& ib = sve.ib;
  int& ic = sve.ic;
  int& j = sve.j;
  int& k = sve.k;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube
  //C       and update the afftected particles through chcell
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  sve.itest = 0;
  //C
  FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ia = i;
        ib = j;
        ic = k;
        if (k >= 1 && k <= 10) {
          if (i == 0) {
            ia = 10;
          }
          if (i == 11) {
            ia = 1;
          }
          if (j == 0) {
            ib = 10;
          }
          if (j == 11) {
            ib = 1;
          }
          chcell(cmn, l, ia, ib, ic, last0, t, tmin, nc);
        }
      }
    }
  }
  //C
}

struct chin3_save
{
  int i;
  int ia;
  int ib;
  int ic;
  int itest;
  int j;
  int k;

  chin3_save() :
    i(fem::int0),
    ia(fem::int0),
    ib(fem::int0),
    ic(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
chin3(
  common& cmn,
  int const& l,
  int const& i1,
  int const& i2,
  int const& i3,
  int const& last0,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chin3);
  // SAVE
  int& i = sve.i;
  int& ia = sve.ia;
  int& ib = sve.ib;
  int& ic = sve.ic;
  int& j = sve.j;
  int& k = sve.k;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube
  //C       and update the afftected particles through chcell
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  sve.itest = 0;
  //C
  FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i == 0) {
          ia = 10;
        }
        else if (i == 11) {
          ia = 1;
        }
        else {
          ia = i;
        }
        if (j == 0) {
          ib = 10;
        }
        else if (j == 11) {
          ib = 1;
        }
        else {
          ib = j;
        }
        if (k == 0) {
          ic = 10;
        }
        else if (k == 11) {
          ic = 1;
        }
        else {
          ic = k;
        }
        chcell(cmn, l, ia, ib, ic, last0, t, tmin, nc);
      }
    }
  }
  //C
}

struct reor_save
{
  int i1;
  int i2;
  int i3;
  int icels0;
  int nc;
  double tmin1;

  reor_save() :
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    icels0(fem::int0),
    nc(fem::int0),
    tmin1(fem::double0)
  {}
};

void
reor(
  common& cmn,
  double const& t,
  double& tmin,
  int const& j,
  int const& last0)
{
  FEM_CMN_SVE(reor);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON ilist1
  const int maxptn = 400001;
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  //
  // SAVE
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& icels0 = sve.icels0;
  int& nc = sve.nc;
  double& tmin1 = sve.tmin1;
  //
  //C       this subroutine is used to fix ct(i) when tm is greater than ct(i)
  //C       next(i) is last1 or last2
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /ilist1/
  //Cd        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
  //Cc      SAVE /ilist5/
  //C
  icels0 = icels(j);
  //C
  i1 = icels0 / 10000;
  i2 = (icels0 - i1 * 10000) / 100;
  i3 = icels0 - i1 * 10000 - i2 * 100;
  //C
  wallc(cmn, j, i1, i2, i3, t, tmin1);
  //C
  if (tmin <= tmin1) {
    nc = last0;
  }
  else {
    tmin = tmin1;
    nc = 0;
  }
  //C
  if (iconfg == 3 || iconfg == 5) {
    chcell(cmn, j, i1, i2, i3, last0, t, tmin, nc);
  }
  else {
    if (i1 == 11 && i2 == 11 && i3 == 11) {
      chout(cmn, j, last0, t, tmin, nc);
    }
    else {
      if (iconfg == 1) {
        chin1(cmn, j, i1, i2, i3, last0, t, tmin, nc);
      }
      else if (iconfg == 2) {
        chin2(cmn, j, i1, i2, i3, last0, t, tmin, nc);
      }
      else if (iconfg == 4) {
        chin3(cmn, j, i1, i2, i3, last0, t, tmin, nc);
      }
    }
  }
  //C
  fixtim(cmn, j, t, tmin1, tmin, nc);
  //C
}

struct dchcel_save
{
  int last0;
  int m;
  int n;
  double tm;

  dchcel_save() :
    last0(fem::int0),
    m(fem::int0),
    n(fem::int0),
    tm(fem::double0)
  {}
};

void
dchcel(
  common& cmn,
  int const& l,
  int const& i,
  int const& j,
  int const& k,
  double const& t)
{
  FEM_CMN_SVE(dchcel);
  // COMMON ilist1
  const int maxptn = 400001;
  arr_cref<int> next(cmn.next, dimension(maxptn));
  arr_cref<int> nic(cmn.nic, dimension(maxptn));
  // COMMON ilist2
  arr_cref<int, 3> icel(cmn.icel, dimension(10, 10, 10));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  //
  // SAVE
  int& last0 = sve.last0;
  int& m = sve.m;
  int& n = sve.n;
  double& tm = sve.tm;
  //
  //C       this subroutine is used to recalculate next collision time for
  //C       particles in the cell i,j,k if the next collision partener is l
  //C
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist5/
  //C
  if (i == 11 || j == 11 || k == 11) {
    if (!(i == 11 && j == 11 && k == 11)) {
      FEM_STOP("cerr");
    }
    m = cmn.icell;
  }
  else {
    m = icel(i, j, k);
  }
  //C
  if (m == 0) {
    return;
  }
  if (next(m) == l) {
    tm = tlarge;
    last0 = 0;
    reor(cmn, t, tm, m, last0);
  }
  n = nic(m);
  if (n == 0) {
    return;
  }
  while (n != m) {
    if (next(n) == l) {
      tm = tlarge;
      last0 = 0;
      reor(cmn, t, tm, n, last0);
    }
    n = nic(n);
  }
  //C
}

struct dchout_save
{
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  int k;
  double td;
  double tt;
  double x1;
  double x2;
  double x3;

  dchout_save() :
    i(fem::int0),
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    j(fem::int0),
    k(fem::int0),
    td(fem::double0),
    tt(fem::double0),
    x1(fem::double0),
    x2(fem::double0),
    x3(fem::double0)
  {}
};

void
dchout(
  common& cmn,
  int const& l,
  int const& ii,
  double const& t)
{
  FEM_CMN_SVE(dchout);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  double& v1 = cmn.v1;
  double& v2 = cmn.v2;
  double& v3 = cmn.v3;
  //
  // SAVE
  int& i = sve.i;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& j = sve.j;
  int& k = sve.k;
  double& td = sve.td;
  double& tt = sve.tt;
  double& x1 = sve.x1;
  double& x2 = sve.x2;
  double& x3 = sve.x3;
  //
  //C       this subroutine is used to check collisions of l with particles when
  //C       l is outside the cube and the collision just happened is a collision
  //C       including a collision with wall (hence we need to use dcheck to throw
  //C       away old collisions that are not in the new neighboring cells.
  //C
  //C       input l,t
  //C
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  tt = ft(l);
  td = t - cmn.size;
  x1 = gx(l) + vx(l) * (t - tt);
  x2 = gy(l) + vy(l) * (t - tt);
  x3 = gz(l) + vz(l) * (t - tt);
  if (td <= 0e0) {
    i1 = integ(x1 / size1) + 6;
    i2 = integ(x2 / size2) + 6;
    i3 = integ(x3 / size3) + 6;
    if (integ(x1 / size1) == x1 / size1 && vx(l) < 0e0) {
      i1 = i1 - 1;
    }
    if (integ(x2 / size2) == x2 / size2 && vy(l) < 0e0) {
      i2 = i2 - 1;
    }
    if (integ(x3 / size3) == x3 / size3 && vz(l) < 0e0) {
      i3 = i3 - 1;
    }
  }
  else {
    i1 = integ(x1 / (size1 + v1 * td)) + 6;
    i2 = integ(x2 / (size2 + v2 * td)) + 6;
    i3 = integ(x3 / (size3 + v3 * td)) + 6;
    //C     10/24/02 (i) below should be (l):
    if (integ(x1 / (size1 + v1 * td)) == x1 / (size1 + v1 * td) && vx(
        l) < (i1 - 6) * v1) {
      i1 = i1 - 1;
    }
    //C     &        vx(i) .lt. (i1 - 6) * v1) i1 = i1 - 1
    if (integ(x2 / (size2 + v2 * td)) == x2 / (size2 + v2 * td) && vy(
        l) < (i2 - 6) * v2) {
      i2 = i2 - 1;
    }
    //C     &        vy(i) .lt. (i2 - 6) * v2) i2 = i2 - 1
    if (integ(x3 / (size3 + v3 * td)) == x3 / (size3 + v3 * td) && vz(
        l) < (i3 - 6) * v3) {
      i3 = i3 - 1;
    }
    //C     &        vz(i) .lt. (i3 - 6) * v3) i3 = i3 - 1
  }
  //C
  if (ii == 1) {
    i = 9;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (j >= 1 && j <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 2) {
    i = 2;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (j >= 1 && j <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 3) {
    j = 9;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i >= 1 && i <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 4) {
    j = 2;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i >= 1 && i <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 5) {
    k = 9;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        if (i >= 1 && i <= 10 && j >= 1 && j <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 6) {
    k = 2;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        if (i >= 1 && i <= 10 && j >= 1 && j <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
}

struct dchin1_save
{
  int i;
  int itest;
  int j;
  int k;

  dchin1_save() :
    i(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
dchin1(
  common& cmn,
  int const& l,
  int const& ii,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t)
{
  FEM_CMN_SVE(dchin1);
  int& i = sve.i;
  int& j = sve.j;
  int& k = sve.k;
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube when the collision just happened is a collision including
  //C       collision with wall
  //C       and update the afftected particles through chkcel
  //C
  //C       input l,ii(specifying the direction of the wall collision),
  //C          i1,i2,i3, (specifying the position of the cell
  //C                    we are going to check)
  //C          t
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  sve.itest = 0;
  //C
  if (ii == 1) {
    if (i1 == 1) {
      goto statement_100;
    }
    if (i1 == 2) {
      if (i2 >= 2 && i2 <= 9 && i3 >= 2 && i3 <= 9) {
        i = 11;
        j = 11;
        k = 11;
        dchcel(cmn, l, i, j, k, t);
      }
      goto statement_100;
    }
    i = i1 - 2;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (j >= 1 && j <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 2) {
    if (i1 == 10) {
      goto statement_100;
    }
    if (i1 == 9) {
      if (i2 >= 2 && i2 <= 9 && i3 >= 2 && i3 <= 9) {
        i = 11;
        j = 11;
        k = 11;
        dchcel(cmn, l, i, j, k, t);
      }
      goto statement_100;
    }
    i = i1 + 2;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (j >= 1 && j <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 3) {
    if (i2 == 1) {
      goto statement_100;
    }
    if (i2 == 2) {
      if (i1 >= 2 && i1 <= 9 && i3 >= 2 && i3 <= 9) {
        i = 11;
        j = 11;
        k = 11;
        dchcel(cmn, l, i, j, k, t);
      }
      goto statement_100;
    }
    j = i2 - 2;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i >= 1 && i <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 4) {
    if (i2 == 10) {
      goto statement_100;
    }
    if (i2 == 9) {
      if (i1 >= 2 && i1 <= 9 && i3 >= 2 && i3 <= 9) {
        i = 11;
        j = 11;
        k = 11;
        dchcel(cmn, l, i, j, k, t);
      }
      goto statement_100;
    }
    j = i2 + 2;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i >= 1 && i <= 10 && k >= 1 && k <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 5) {
    if (i3 == 1) {
      goto statement_100;
    }
    if (i3 == 2) {
      if (i1 >= 2 && i1 <= 9 && i2 >= 2 && i2 <= 9) {
        i = 11;
        j = 11;
        k = 11;
        dchcel(cmn, l, i, j, k, t);
      }
      goto statement_100;
    }
    k = i3 - 2;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        if (i >= 1 && i <= 10 && j >= 1 && j <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  if (ii == 6) {
    if (i3 == 10) {
      goto statement_100;
    }
    if (i3 == 9) {
      if (i1 >= 2 && i1 <= 9 && i2 >= 2 && i2 <= 9) {
        i = 11;
        j = 11;
        k = 11;
        dchcel(cmn, l, i, j, k, t);
      }
      goto statement_100;
    }
    k = i3 + 2;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        if (i >= 1 && i <= 10 && j >= 1 && j <= 10) {
          dchcel(cmn, l, i, j, k, t);
        }
      }
    }
  }
  //C
  statement_100:;
  //C
}

struct dchin2_save
{
  int i;
  int ia;
  int ib;
  int ic;
  int j;
  int k;

  dchin2_save() :
    i(fem::int0),
    ia(fem::int0),
    ib(fem::int0),
    ic(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
dchin2(
  common& cmn,
  int const& l,
  int const& ii,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t)
{
  FEM_CMN_SVE(dchin2);
  int& i = sve.i;
  int& ia = sve.ia;
  int& ib = sve.ib;
  int& ic = sve.ic;
  int& j = sve.j;
  int& k = sve.k;
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube when the collision just happened is a collision including
  //C       collision with wall
  //C       and update the afftected particles through chkcel
  //C
  //C       input l,ii(specifying the direction of the wall collision),
  //C          i1,i2,i3, (specifying the position of the cell
  //C                    we are going to check)
  //C          t
  //C
  if (ii == 1) {
    i = i1 - 2;
    if (i <= 0) {
      i += 10;
    }
    ia = i;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ib = j;
        ic = k;
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        if (k >= 1 && k <= 10) {
          dchcel(cmn, l, ia, ib, ic, t);
        }
      }
    }
  }
  //C
  if (ii == 2) {
    i = i1 + 2;
    if (i >= 11) {
      i = i - 10;
    }
    ia = i;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ib = j;
        ic = k;
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        if (k >= 1 && k <= 10) {
          dchcel(cmn, l, ia, ib, ic, t);
        }
      }
    }
  }
  //C
  if (ii == 3) {
    j = i2 - 2;
    if (j <= 0) {
      j += 10;
    }
    ib = j;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ia = i;
        ic = k;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (k >= 1 && k <= 10) {
          dchcel(cmn, l, ia, ib, ic, t);
        }
      }
    }
  }
  //C
  if (ii == 4) {
    j = i2 + 2;
    if (j >= 11) {
      j = j - 10;
    }
    ib = j;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ia = i;
        ic = k;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (k >= 1 && k <= 10) {
          dchcel(cmn, l, ia, ib, ic, t);
        }
      }
    }
  }
  //C
  if (ii == 5) {
    if (i3 == 2) {
      goto statement_100;
    }
    k = i3 - 2;
    ic = k;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        ia = i;
        ib = j;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  if (ii == 6) {
    if (i3 == 9) {
      goto statement_100;
    }
    k = i3 + 2;
    ic = k;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        ia = i;
        ib = j;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  statement_100:;
  //C
}

struct dchin3_save
{
  int i;
  int ia;
  int ib;
  int ic;
  int j;
  int k;

  dchin3_save() :
    i(fem::int0),
    ia(fem::int0),
    ib(fem::int0),
    ic(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
dchin3(
  common& cmn,
  int const& l,
  int const& ii,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t)
{
  FEM_CMN_SVE(dchin3);
  // SAVE
  int& i = sve.i;
  int& ia = sve.ia;
  int& ib = sve.ib;
  int& ic = sve.ic;
  int& j = sve.j;
  int& k = sve.k;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube when the collision just happened is a collision including
  //C       collision with wall
  //C       and update the afftected particles through chkcel
  //C
  //C       input l,ii(specifying the direction of the wall collision),
  //C          i1,i2,i3, (specifying the position of the cell
  //C                    we are going to check)
  //C          t
  //C
  if (ii == 1) {
    i = i1 - 2;
    if (i <= 0) {
      i += 10;
    }
    ia = i;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ib = j;
        ic = k;
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        if (k == 0) {
          ic = 10;
        }
        if (k == 11) {
          ic = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  if (ii == 2) {
    i = i1 + 2;
    if (i >= 11) {
      i = i - 10;
    }
    ia = i;
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ib = j;
        ic = k;
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        if (k == 0) {
          ic = 10;
        }
        if (k == 11) {
          ic = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  if (ii == 3) {
    j = i2 - 2;
    if (j <= 0) {
      j += 10;
    }
    ib = j;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ia = i;
        ic = k;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (k == 0) {
          ic = 10;
        }
        if (k == 11) {
          ic = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  if (ii == 4) {
    j = i2 + 2;
    if (j >= 11) {
      j = j - 10;
    }
    ib = j;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ia = i;
        ic = k;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (k == 0) {
          ic = 10;
        }
        if (k == 11) {
          ic = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  if (ii == 5) {
    k = i3 - 2;
    if (k <= 0) {
      k += 10;
    }
    ic = k;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        ia = i;
        ib = j;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
  if (ii == 6) {
    k = i3 + 2;
    if (k >= 11) {
      k = k - 10;
    }
    ic = k;
    FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
      FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
        ia = i;
        ib = j;
        if (i == 0) {
          ia = 10;
        }
        if (i == 11) {
          ia = 1;
        }
        if (j == 0) {
          ib = 10;
        }
        if (j == 11) {
          ib = 1;
        }
        dchcel(cmn, l, ia, ib, ic, t);
      }
    }
  }
  //C
}

struct cellre_save
{
  double ctmp;
  double ddt;
  double dtt;
  bool good;
  int i1;
  int i2;
  int i3;
  int icels0;
  int ii;
  int j;
  int k;
  double otmp;
  double t0;
  double tmin1;

  cellre_save() :
    ctmp(fem::double0),
    ddt(fem::double0),
    dtt(fem::double0),
    good(fem::bool0),
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    icels0(fem::int0),
    ii(fem::int0),
    j(fem::int0),
    k(fem::int0),
    otmp(fem::double0),
    t0(fem::double0),
    tmin1(fem::double0)
  {}
};

void
cellre(
  common& cmn,
  int const& i,
  double const& t)
{
  FEM_CMN_SVE(cellre);
  int& iconfg = cmn.iconfg;
  const int maxptn = 400001;
  arr_ref<double> gx(cmn.gx, dimension(maxptn));
  arr_ref<double> gy(cmn.gy, dimension(maxptn));
  arr_ref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  arr_ref<double> dgxa(cmn.dgxa, dimension(maxptn));
  arr_ref<double> dgya(cmn.dgya, dimension(maxptn));
  arr_ref<double> dgza(cmn.dgza, dimension(maxptn));
  arr_ref<int> next(cmn.next, dimension(maxptn));
  arr_ref<int> icsta(cmn.icsta, dimension(maxptn));
  arr_cref<int> nic(cmn.nic, dimension(maxptn));
  arr_ref<int> icels(cmn.icels, dimension(maxptn));
  int& icell = cmn.icell;
  arr_ref<int, 3> icel(cmn.icel, dimension(10, 10, 10));
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  int& ichkpt = cmn.ichkpt;
  arr_ref<double> ct(cmn.ct, dimension(maxptn));
  arr_ref<double> ot(cmn.ot, dimension(maxptn));
  //
  double& ctmp = sve.ctmp;
  double& ddt = sve.ddt;
  double& dtt = sve.dtt;
  bool& good = sve.good;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& icels0 = sve.icels0;
  int& ii = sve.ii;
  int& j = sve.j;
  int& k = sve.k;
  double& otmp = sve.otmp;
  double& t0 = sve.t0;
  double& tmin1 = sve.tmin1;
  //C       this subroutine is used for changing the cell of a particle that
  //C       collide with the wall
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //C
  //C       this happens before update the /prec2/ common; in contrast with
  //C       scat which happens after updating the glue common
  //C
  t0 = t;
  //C
  statement_1000:
  //C
  if (iconfg == 3 || iconfg == 5) {
    k = fem::mod(icsta(i), 10);
    //C
    if (k == 1) {
      gx(i) = gx(i) - 10e0 * size1;
      dgxa(i) += 10e0 * size1;
      FEM_DO_SAFE(ii, 1, ichkpt) {
        if (next(ii) == i) {
          dgxa(ii) = dgxa(ii) - 10e0 * size1;
        }
      }
    }
    if (k == 2) {
      gx(i) += 10e0 * size1;
      dgxa(i) = dgxa(i) - 10e0 * size1;
      FEM_DO_SAFE(ii, 1, ichkpt) {
        if (next(ii) == i) {
          dgxa(ii) += 10e0 * size1;
        }
      }
    }
    if (k == 3) {
      gy(i) = gy(i) - 10e0 * size2;
      dgya(i) += 10e0 * size2;
      FEM_DO_SAFE(ii, 1, ichkpt) {
        if (next(ii) == i) {
          dgya(ii) = dgya(ii) - 10e0 * size2;
        }
      }
    }
    if (k == 4) {
      gy(i) += 10e0 * size2;
      dgya(i) = dgya(i) - 10e0 * size2;
      FEM_DO_SAFE(ii, 1, ichkpt) {
        if (next(ii) == i) {
          dgya(ii) += 10e0 * size2;
        }
      }
    }
    if (iconfg == 5) {
      if (k == 5) {
        gz(i) = gz(i) - 10e0 * size3;
        dgza(i) += 10e0 * size3;
        FEM_DO_SAFE(ii, 1, ichkpt) {
          if (next(ii) == i) {
            dgza(ii) = dgza(ii) - 10e0 * size3;
          }
        }
      }
      if (k == 6) {
        gz(i) += 10e0 * size3;
        dgza(i) = dgza(i) - 10e0 * size3;
        FEM_DO_SAFE(ii, 1, ichkpt) {
          if (next(ii) == i) {
            dgza(ii) += 10e0 * size3;
          }
        }
      }
    }
  }
  else {
    icels0 = icels(i);
    //C
    i1 = icels0 / 10000;
    i2 = (icels0 - i1 * 10000) / 100;
    i3 = icels0 - i1 * 10000 - i2 * 100;
    //C
    //Cc       for particle inside the cube
    if (i1 >= 1 && i1 <= 10 && i2 >= 1 && i2 <= 10 && i3 >= 1 && i3 <= 10) {
      //C
      //C       this assignment takes care of nic(i)=0 automatically
      if (icel(i1, i2, i3) == i) {
        icel(i1, i2, i3) = nic(i);
      }
      //C
      //C1      rearrange the old cell
      //C
      oldcre(cmn, i);
      //C
      //C2      rearrange the new cell
      //C
      k = fem::mod(icsta(i), 10);
      //C
      //C2.1    particle goes out of the cube
      if (iconfg == 1) {
        good = (i1 == 1 && k == 2) || (i1 == 10 && k == 1) || (
          i2 == 1 && k == 4) || (i2 == 10 && k == 3) || (i3 == 1 &&
          k == 6) || (i3 == 10 && k == 5);
      }
      if (iconfg == 2) {
        good = (i3 == 1 && k == 6) || (i3 == 10 && k == 5);
      }
      if (good) {
        //C
        //C                j = icell
        newcre(cmn, i, icell);
        //C                 icell = j
        //C
        icels(i) = 111111;
        //C
        //C2.2    particle moves inside the cube
      }
      else {
        //C
        if (k == 1) {
          i1++;
        }
        if (k == 2) {
          i1 = i1 - 1;
        }
        if (k == 3) {
          i2++;
        }
        if (k == 4) {
          i2 = i2 - 1;
        }
        if (k == 5) {
          i3++;
        }
        if (k == 6) {
          i3 = i3 - 1;
        }
        //C
        if (iconfg == 2 || iconfg == 4) {
          if (i1 == 0) {
            i1 = 10;
            gx(i) += 10e0 * size1;
          }
          if (i1 == 11) {
            i1 = 1;
            gx(i) = gx(i) - 10e0 * size1;
          }
          if (i2 == 0) {
            i2 = 10;
            gy(i) += 10e0 * size2;
          }
          if (i2 == 11) {
            i2 = 1;
            gy(i) = gy(i) - 10e0 * size2;
          }
          if (iconfg == 4) {
            if (i3 == 0) {
              i3 = 10;
              gz(i) += 10e0 * size3;
            }
            if (i3 == 11) {
              i3 = 1;
              gz(i) = gz(i) - 10e0 * size3;
            }
          }
        }
        //C
        j = icel(i1, i2, i3);
        //C
        newcre(cmn, i, j);
        //C       in case icel changes
        //C
        icel(i1, i2, i3) = j;
        //C
        icels(i) = i1 * 10000 + i2 * 100 + i3;
        //C
      }
      //C
      //Cc       for particles outside the cube
    }
    else {
      //C
      if (icell == i) {
        icell = nic(i);
      }
      //C
      oldcre(cmn, i);
      //C
      k = fem::mod(icsta(i), 10);
      //C
      ddt = t - ft(i);
      dtt = t - cmn.size;
      if (dtt <= 0e0) {
        i1 = integ((gx(i) + vx(i) * ddt) / size1) + 6;
        i2 = integ((gy(i) + vy(i) * ddt) / size2) + 6;
        i3 = integ((gz(i) + vz(i) * ddt) / size3) + 6;
      }
      else {
        i1 = integ((gx(i) + vx(i) * ddt) / (size1 + cmn.v1 * dtt)) + 6;
        i2 = integ((gy(i) + vy(i) * ddt) / (size2 + cmn.v2 * dtt)) + 6;
        i3 = integ((gz(i) + vz(i) * ddt) / (size3 + cmn.v3 * dtt)) + 6;
      }
      //C
      if (k == 1) {
        i1 = 1;
      }
      if (k == 2) {
        i1 = 10;
      }
      if (k == 3) {
        i2 = 1;
      }
      if (k == 4) {
        i2 = 10;
      }
      if (k == 5) {
        i3 = 1;
      }
      if (k == 6) {
        i3 = 10;
      }
      //C
      j = icel(i1, i2, i3);
      newcre(cmn, i, j);
      icel(i1, i2, i3) = j;
      //C
      icels(i) = i1 * 10000 + i2 * 100 + i3;
      //C
    }
  }
  //C
  if (next(i) != 0) {
    otmp = ot(next(i));
    ctmp = ct(next(i));
  }
  //C
  if (i1 == 11 && i2 == 11 && i3 == 11) {
    dchout(cmn, i, k, t);
  }
  else {
    if (iconfg == 1) {
      dchin1(cmn, i, k, i1, i2, i3, t);
    }
    else if (iconfg == 2) {
      dchin2(cmn, i, k, i1, i2, i3, t);
    }
    else if (iconfg == 4) {
      dchin3(cmn, i, k, i1, i2, i3, t);
    }
  }
  //C
  if (icsta(i) / 10 == 11) {
    ot(next(i)) = otmp;
    ct(next(i)) = ctmp;
    next(next(i)) = i;
    wallc(cmn, i, i1, i2, i3, t0, tmin1);
    if (tmin1 < ct(i)) {
      icsta(i) += 10;
      t0 = tmin1;
      goto statement_1000;
    }
  }
  //C
}

struct newpos_save
{
  double dt1;

  newpos_save() :
    dt1(fem::double0)
  {}
};

void
newpos(
  common& cmn,
  double const& /* t */,
  int const& i)
{
  FEM_CMN_SVE(newpos);
  // COMMON prec2
  const int maxptn = 400001;
  arr_ref<double> gx(cmn.gx, dimension(maxptn));
  arr_ref<double> gy(cmn.gy, dimension(maxptn));
  arr_ref<double> gz(cmn.gz, dimension(maxptn));
  arr_ref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec4
  arr_cref<double> vx(cmn.vx, dimension(maxptn));
  arr_cref<double> vy(cmn.vy, dimension(maxptn));
  arr_cref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON prec5
  arr_ref<double> eta(cmn.eta, dimension(maxptn));
  arr_ref<double> tau(cmn.tau, dimension(maxptn));
  // COMMON ilist5
  arr_cref<double> ct(cmn.ct, dimension(maxptn));
  //
  // SAVE
  double& dt1 = sve.dt1;
  //
  //C
  //C       this subroutine is used to calculate the 2 particle scattering
  //C       get new position
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /prec5/
  //Cc      SAVE /ilist5/
  //C
  dt1 = ct(i) - ft(i);
  //C
  gx(i) += vx(i) * dt1;
  gy(i) += vy(i) * dt1;
  gz(i) += vz(i) * dt1;
  ft(i) = ct(i);
  //C
  if (cmn.iconfg <= 3) {
    if (ft(i) <= fem::abs(gz(i))) {
      eta(i) = 1000000.e0;
    }
    else {
      eta(i) = 0.5e0 * fem::log((ft(i) + gz(i)) / (ft(i) - gz(i)));
    }
    //Clin-8/2015 to avoid IEEE_OVERFLOW_FLAG:
    //C           tau(i) = ft(i) / cosh(eta(i))
    if (eta(i) < 1000000.e0) {
      tau(i) = ft(i) / fem::cosh(eta(i));
    }
    else {
      tau(i) = 1e-10;
    }
    //C
  }
  //C
}

struct getht_save
{
  int iseed;
  double rx;
  double xm2;
  double xmp2;
  double xmu2;

  getht_save() :
    iseed(fem::int0),
    rx(fem::double0),
    xm2(fem::double0),
    xmp2(fem::double0),
    xmu2(fem::double0)
  {}
};

void
getht(
  common& cmn,
  int const& /* iscat */,
  int const& /* jscat */,
  double const& pp2,
  double& that)
{
  FEM_CMN_SVE(getht);
  // SAVE
  int& iseed = sve.iseed;
  double& rx = sve.rx;
  double& xm2 = sve.xm2;
  double& xmp2 = sve.xmp2;
  double& xmu2 = sve.xmu2;
  //
  //C
  //C       this subroutine is used to get \hat{t} for a particular processes
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /anim/
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  const double hbarc = 0.197327054e0;
  xmu2 = fem::pow2((hbarc * cmn.xmu));
  xmp2 = fem::pow2(cmn.xmp);
  xm2 = xmu2 + xmp2;
  rx = ran1(cmn, iseed);
  that = xm2 * (1e0 + 1e0 / ((1e0 - xm2 / (4e0 * pp2 + xm2)) * rx - 1e0));
  //Ctest off isotropic scattering:
  //C     &     + 1d0/((1d0 - xm2 / (4d0 * pp2 + xm2)) * ran1(2) - 1d0))
  //C        if(izpc.eq.100) that=-4d0*pp2*ran1(2)
  if (cmn.izpc == 100) {
    that = -4e0 * pp2 * rx;
  }
  //C
}

struct cropro_save
{
  fem::variant_bindings cprod_bindings;
};

void
cropro(
  common& cmn,
  double const& vx1,
  double const& vy1,
  double const& vz1,
  double const& vx2,
  double const& vy2,
  double const& vz2)
{
  FEM_CMN_SVE(cropro);
  common_variant cprod(cmn.common_cprod, sve.cprod_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<double> vx3;
      mbr<double> vy3;
      mbr<double> vz3;
      cprod.allocate(), vx3, vy3, vz3;
    }
  }
  double& vx3 = cprod.bind<double>();
  double& vy3 = cprod.bind<double>();
  double& vz3 = cprod.bind<double>();
  //C
  //C     this subroutine is used to calculate the cross product of
  //C     (vx1,vy1,vz1) and (vx2,vy2,vz2) and get the result (vx3,vy3,vz3)
  //C     and put the vector into common /cprod/
  //C
  //Cc      SAVE /cprod/
  //C
  vx3 = vy1 * vz2 - vz1 * vy2;
  vy3 = vz1 * vx2 - vx1 * vz2;
  vz3 = vx1 * vy2 - vy1 * vx2;
  //C
}

struct xnormv_save
{
  double vv;

  xnormv_save() :
    vv(fem::double0)
  {}
};

void
xnormv(
  common& cmn,
  double& vx,
  double& vy,
  double& vz)
{
  FEM_CMN_SVE(xnormv);
  // SAVE
  double& vv = sve.vv;
  //
  //C
  //C      this subroutine is used to get a normalized vector
  //C
  //Clin-7/20/01:
  //C      vv = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
  vv = fem::dsqrt(fem::pow2(vx) + fem::pow2(vy) + fem::pow2(vz));
  vx = vx / vv;
  vy = vy / vv;
  vz = vz / vv;
  //C
}

struct zprota_save
{
  double a11;
  double a12;
  double a13;
  double a21;
  double a22;
  double a23;
  double a31;
  double a32;
  double a33;
  double c;
  double omc;
  double s;
  double vx;
  double vy;
  double vz;

  zprota_save() :
    a11(fem::double0),
    a12(fem::double0),
    a13(fem::double0),
    a21(fem::double0),
    a22(fem::double0),
    a23(fem::double0),
    a31(fem::double0),
    a32(fem::double0),
    a33(fem::double0),
    c(fem::double0),
    omc(fem::double0),
    s(fem::double0),
    vx(fem::double0),
    vy(fem::double0),
    vz(fem::double0)
  {}
};

//C
//Cbz1/29/99
//C      subroutine rotate(xn1, xn2, xn3, theta, v1, v2, v3)
void
zprota(
  common& cmn,
  double const& xn1,
  double const& xn2,
  double const& xn3,
  double const& theta,
  double& v1,
  double& v2,
  double& v3)
{
  FEM_CMN_SVE(zprota);
  // SAVE
  double& a11 = sve.a11;
  double& a12 = sve.a12;
  double& a13 = sve.a13;
  double& a21 = sve.a21;
  double& a22 = sve.a22;
  double& a23 = sve.a23;
  double& a31 = sve.a31;
  double& a32 = sve.a32;
  double& a33 = sve.a33;
  double& c = sve.c;
  double& omc = sve.omc;
  double& s = sve.s;
  double& vx = sve.vx;
  double& vy = sve.vy;
  double& vz = sve.vz;
  //
  //Cbz1/29/99end
  //C
  //C     this subroutine is used to rotate the vector (v1,v2,v3) by an angle theta
  //C     around the unit vector (xn1, xn2, xn3)
  //C
  vx = v1;
  vy = v2;
  vz = v3;
  c = fem::cos(theta);
  omc = 1e0 - c;
  s = fem::sin(theta);
  a11 = fem::pow2(xn1) * omc + c;
  a12 = xn1 * xn2 * omc - s * xn3;
  a13 = xn1 * xn3 * omc + s * xn2;
  a21 = xn1 * xn2 * omc + s * xn3;
  a22 = fem::pow2(xn2) * omc + c;
  a23 = xn2 * xn3 * omc - s * xn1;
  a31 = xn1 * xn3 * omc - s * xn2;
  a32 = xn3 * xn2 * omc + s * xn1;
  a33 = fem::pow2(xn3) * omc + c;
  v1 = vx * a11 + vy * a12 + vz * a13;
  v2 = vx * a21 + vy * a22 + vz * a23;
  v3 = vx * a31 + vy * a32 + vz * a33;
  //C
}

struct newmom_save
{
  fem::variant_bindings cprod_bindings;
  double bex;
  double bey;
  double bez;
  double e1;
  double e2;
  int i1;
  int i2;
  int icels1;
  int icels2;
  int j1;
  int j2;
  int k1;
  int k2;
  double pp2;
  double px1;
  double px2;
  double py1;
  double py2;
  double pz1;
  double pz2;
  double rap1;
  double rap2;
  double rts2;
  double t1;
  double t2;
  double that;
  double theta;
  double x1;
  double x2;
  double y1;
  double y2;
  double z1;
  double z2;

  newmom_save() :
    bex(fem::double0),
    bey(fem::double0),
    bez(fem::double0),
    e1(fem::double0),
    e2(fem::double0),
    i1(fem::int0),
    i2(fem::int0),
    icels1(fem::int0),
    icels2(fem::int0),
    j1(fem::int0),
    j2(fem::int0),
    k1(fem::int0),
    k2(fem::int0),
    pp2(fem::double0),
    px1(fem::double0),
    px2(fem::double0),
    py1(fem::double0),
    py2(fem::double0),
    pz1(fem::double0),
    pz2(fem::double0),
    rap1(fem::double0),
    rap2(fem::double0),
    rts2(fem::double0),
    t1(fem::double0),
    t2(fem::double0),
    that(fem::double0),
    theta(fem::double0),
    x1(fem::double0),
    x2(fem::double0),
    y1(fem::double0),
    y2(fem::double0),
    z1(fem::double0),
    z2(fem::double0)
  {}
};

void
newmom(
  common& cmn,
  double const& /* t */)
{
  FEM_CMN_SVE(newmom);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON para6
  double& centy = cmn.centy;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_ref<double> px(cmn.px, dimension(maxptn));
  arr_ref<double> py(cmn.py, dimension(maxptn));
  arr_ref<double> pz(cmn.pz, dimension(maxptn));
  arr_ref<double> e(cmn.e, dimension(maxptn));
  arr_cref<double> xmass(cmn.xmass, dimension(maxptn));
  // COMMON prec4
  arr_ref<double> vx(cmn.vx, dimension(maxptn));
  arr_ref<double> vy(cmn.vy, dimension(maxptn));
  arr_ref<double> vz(cmn.vz, dimension(maxptn));
  // COMMON prec5
  arr_ref<double> rap(cmn.rap, dimension(maxptn));
  // COMMON aurec2
  arr_cref<double> dgxa(cmn.dgxa, dimension(maxptn));
  arr_cref<double> dgya(cmn.dgya, dimension(maxptn));
  arr_cref<double> dgza(cmn.dgza, dimension(maxptn));
  // COMMON ilist1
  int& iscat = cmn.iscat;
  int& jscat = cmn.jscat;
  arr_ref<int> last(cmn.last, dimension(maxptn));
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  // COMMON lor
  double& enenew = cmn.enenew;
  double& pxnew = cmn.pxnew;
  double& pynew = cmn.pynew;
  double& pznew = cmn.pznew;
  // COMMON rndm2
  int& iff = cmn.iff;
  // COMMON frzprc
  arr_cref<int> ifrz(cmn.ifrz, dimension(maxptn));
  //
  common_variant cprod(cmn.common_cprod, sve.cprod_bindings);
  // SAVE
  double& bex = sve.bex;
  double& bey = sve.bey;
  double& bez = sve.bez;
  double& e1 = sve.e1;
  double& e2 = sve.e2;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& icels1 = sve.icels1;
  int& icels2 = sve.icels2;
  int& j1 = sve.j1;
  int& j2 = sve.j2;
  int& k1 = sve.k1;
  int& k2 = sve.k2;
  double& pp2 = sve.pp2;
  double& px1 = sve.px1;
  double& px2 = sve.px2;
  double& py1 = sve.py1;
  double& py2 = sve.py2;
  double& pz1 = sve.pz1;
  double& pz2 = sve.pz2;
  double& rap1 = sve.rap1;
  double& rap2 = sve.rap2;
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  double& that = sve.that;
  double& theta = sve.theta;
  double& x1 = sve.x1;
  double& x2 = sve.x2;
  double& y1 = sve.y1;
  double& y2 = sve.y2;
  double& z1 = sve.z1;
  double& z2 = sve.z2;
  //
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<double> xn1;
      mbr<double> xn2;
      mbr<double> xn3;
      cprod.allocate(), xn1, xn2, xn3;
    }
  }
  double& xn1 = cprod.bind<double>();
  double& xn2 = cprod.bind<double>();
  double& xn3 = cprod.bind<double>();
  //C
  //C       this subroutine is used to calculate the 2 particle scattering
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Ctrans
  //Cc      SAVE /para6/
  //Ctransend
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /prec5/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /lor/
  //Cc      SAVE /cprod/
  //Cc      SAVE /rndm2/
  //Cc      SAVE /anim/
  //Cc      SAVE /frzprc/
  //C
  //Clin-6/06/02 no momentum change for partons already frozen out,
  //C     however, spatial upgrade is needed to ensure overall system freezeout:
  if (cmn.isoft == 5) {
    if (ifrz(iscat) == 1 || ifrz(jscat) == 1) {
      last(iscat) = jscat;
      last(jscat) = iscat;
      return;
    }
  }
  //Clin-6/06/02-end
  //C
  //C       iff is used to randomize the interaction to have both attractive and
  //C        repulsive
  //C
  iff = -iff;
  //C
  if (iconfg == 2 || iconfg == 4) {
    icels1 = icels(iscat);
    i1 = icels1 / 10000;
    j1 = (icels1 - i1 * 10000) / 100;
    icels2 = icels(jscat);
    i2 = icels2 / 10000;
    j2 = (icels2 - i2 * 10000) / 100;
    if (iconfg == 4) {
      k1 = icels1 - i1 * 10000 - j1 * 100;
      k2 = icels2 - i2 * 10000 - j2 * 100;
    }
  }
  //C
  px1 = px(iscat);
  py1 = py(iscat);
  pz1 = pz(iscat);
  e1 = e(iscat);
  x1 = gx(iscat);
  y1 = gy(iscat);
  z1 = gz(iscat);
  t1 = ft(iscat);
  px2 = px(jscat);
  py2 = py(jscat);
  pz2 = pz(jscat);
  e2 = e(jscat);
  //C
  if (iconfg == 1) {
    x2 = gx(jscat);
    y2 = gy(jscat);
    z2 = gz(jscat);
  }
  else if (iconfg == 2 || iconfg == 4) {
    if (i1 - i2 > 5) {
      x2 = gx(jscat) + 10e0 * size1;
    }
    else if (i1 - i2 <  - 5) {
      x2 = gx(jscat) - 10e0 * size1;
    }
    else {
      x2 = gx(jscat);
    }
    if (j1 - j2 > 5) {
      y2 = gy(jscat) + 10e0 * size2;
    }
    else if (j1 - j2 <  - 5) {
      y2 = gy(jscat) - 10e0 * size2;
    }
    else {
      y2 = gy(jscat);
    }
    if (iconfg == 4) {
      if (k1 - k2 > 5) {
        z2 = gz(jscat) + 10e0 * size3;
      }
      else if (k1 - k2 <  - 5) {
        z2 = gz(jscat) - 10e0 * size3;
      }
      else {
        z2 = gz(jscat);
      }
    }
    else {
      z2 = gz(jscat);
    }
  }
  else if (iconfg == 3 || iconfg == 5) {
    x2 = gx(jscat) + dgxa(jscat);
    y2 = gy(jscat) + dgya(jscat);
    if (iconfg == 5) {
      z2 = gz(jscat) + dgza(jscat);
    }
    else {
      z2 = gz(jscat);
    }
  }
  t2 = ft(jscat);
  //Ctrans
  sve.rts2 = fem::pow2((e1 + e2)) - fem::pow2((px1 + px2)) -
    fem::pow2((py1 + py2)) - fem::pow2((pz1 + pz2));
  //Ctransend
  bex = (px1 + px2) / (e1 + e2);
  bey = (py1 + py2) / (e1 + e2);
  bez = (pz1 + pz2) / (e1 + e2);
  //C
  //Clin-11/2015-ctest off
  //C        write(99,*) 'iscat,jscat,etotalA=',iscat,jscat,e1+e2
  //C
  lorenz(cmn, e1, px1, py1, pz1, bex, bey, bez);
  //Cc      SAVE pxnew, ..., values for later use.
  px1 = pxnew;
  py1 = pynew;
  pz1 = pznew;
  e1 = enenew;
  //C
  pp2 = fem::pow2(pxnew) + fem::pow2(pynew) + fem::pow2(pznew);
  getht(cmn, iscat, jscat, pp2, that);
  theta = fem::dacos(that / (2e0 * pp2) + 1e0);
  theta = fem::dble(iff) * theta;
  //C
  //C       we boost to the cm frame, get rotation axis, and rotate 1 particle
  //C       momentum
  //C
  lorenz(cmn, t1, x1, y1, z1, bex, bey, bez);
  //C
  x1 = pxnew;
  y1 = pynew;
  z1 = pznew;
  //C
  lorenz(cmn, t2, x2, y2, z2, bex, bey, bez);
  //C
  x2 = pxnew;
  y2 = pynew;
  z2 = pznew;
  //C
  //C       notice now pxnew, ..., are new positions
  cropro(cmn, x1 - x2, y1 - y2, z1 - z2, px1, py1, pz1);
  //C
  xnormv(cmn, xn1, xn2, xn3);
  //C
  //Cbz1/29/99
  //C        call rotate(xn1, xn2, xn3, theta, px1, py1, pz1)
  zprota(cmn, xn1, xn2, xn3, theta, px1, py1, pz1);
  //Cbz1/29/99end
  //C
  //C       we invert the momentum to get the other particle's momentum
  px2 = -px1;
  py2 = -py1;
  pz2 = -pz1;
  //Clin-4/13/01: modify in case m1, m2 are different:
  //C        e2 = e1
  e2 = fem::dsqrt(fem::pow2(px2) + fem::pow2(py2) + fem::pow2(pz2) +
    fem::pow2(xmass(jscat)));
  //C
  //Clin-11/2015-ctest off
  //C        write(99,*) 'iscat,jscat,masses= ',iscat,jscat,
  //C     1       xmass(iscat),xmass(jscat)
  //C
  //C       boost the 2 particle 4 momentum back to lab frame
  lorenz(cmn, e1, px1, py1, pz1, -bex, -bey, -bez);
  px(iscat) = pxnew;
  py(iscat) = pynew;
  pz(iscat) = pznew;
  e(iscat) = enenew;
  lorenz(cmn, e2, px2, py2, pz2, -bex, -bey, -bez);
  px(jscat) = pxnew;
  py(jscat) = pynew;
  pz(jscat) = pznew;
  e(jscat) = enenew;
  //C
  //Clin-11/2015-ctest off
  //C        write(99,*) 'iscat,jscat,etotalB= ',iscat,jscat,
  //C     1       e(iscat)+e(jscat)
  //C
  vx(iscat) = px(iscat) / e(iscat);
  vy(iscat) = py(iscat) / e(iscat);
  vz(iscat) = pz(iscat) / e(iscat);
  vx(jscat) = px(jscat) / e(jscat);
  vy(jscat) = py(jscat) / e(jscat);
  vz(jscat) = pz(jscat) / e(jscat);
  //C
  last(iscat) = jscat;
  last(jscat) = iscat;
  //C
  if (iconfg <= 3) {
    if (e(iscat) <= fem::abs(pz(iscat))) {
      rap(iscat) = 1000000.e0;
    }
    else {
      rap(iscat) = 0.5e0 * fem::log((e(iscat) + pz(iscat)) / (e(
        iscat) - pz(iscat)));
    }
    //C
    if (e(jscat) <= fem::abs(pz(jscat))) {
      rap(jscat) = 1000000.e0;
    }
    else {
      rap(jscat) = 0.5e0 * fem::log((e(jscat) + pz(jscat)) / (e(
        jscat) - pz(jscat)));
    }
    //C
    //Ctrans
    rap1 = rap(iscat);
    rap2 = rap(jscat);
    //C
    if ((rap1 < centy + 0.5e0 && rap1 > centy - 0.5e0)) {
      //C              write (9, *) sqrt(ft(iscat) ** 2 - gz(iscat) ** 2), rts2
    }
    if ((rap2 < centy + 0.5e0 && rap2 > centy - 0.5e0)) {
      //C              write (9, *) sqrt(ft(jscat) ** 2 - gz(jscat) ** 2), rts2
    }
    //Ctransend
  }
  //C
  //Clin-11/2015-ctest off
  //C        write(99,*) 'iscat,jscat,xmp,xmu,that=',iscat,jscat,xmp,xmu,that
  //C
}

void
scat(
  common& cmn,
  double const& t,
  int const& iscat,
  int const& jscat)
{
  //C
  //C       this subroutine is used to calculate the 2 particle scattering
  //C
  newpos(cmn, t, iscat);
  newpos(cmn, t, jscat);
  newmom(cmn, t);
  //C
}

void
ck(
  common& cmn,
  int const& l,
  int& ick)
{
  // COMMON ilist1
  int& iscat = cmn.iscat;
  int& jscat = cmn.jscat;
  int& ictype = cmn.ictype;
  // COMMON ilist4
  int& ifmpt = cmn.ifmpt;
  //
  //C       this subroutine is used for chcell to check whether l should be
  //C       checked or not for updating tmin, nc
  //C       input l
  //C       output ick
  //C       if ick=1, l should be checked, otherwise it should not be.
  //C
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist4/
  //C
  ick = 1;
  if (ictype == 1) {
    if (l == ifmpt) {
      ick = 0;
    }
  }
  else if (ictype == 0 || ictype == 3 || ictype == 4) {
    if (l == iscat || l == jscat) {
      ick = 0;
    }
  }
  else {
    if (l == iscat || l == jscat || l == ifmpt) {
      ick = 0;
    }
  }
  //C       notice il is either iscat or jscat, or ifmpt, we deal with them
  //C       seperately according to ictype
  //C
}

struct ud2_save
{
  bool allok;
  int i1;
  int i2;
  int i3;
  int icels0;
  double t1;
  double t2;
  double tm;
  double tmin1;

  ud2_save() :
    allok(fem::bool0),
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    icels0(fem::int0),
    t1(fem::double0),
    t2(fem::double0),
    tm(fem::double0),
    tmin1(fem::double0)
  {}
};

void
ud2(
  common& cmn,
  int const& i,
  int const& j,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(ud2);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON aurec1
  int& jxa = cmn.jxa;
  int& jya = cmn.jya;
  int& jza = cmn.jza;
  // COMMON aurec2
  const int maxptn = 400001;
  arr_ref<double> dgxa(cmn.dgxa, dimension(maxptn));
  arr_ref<double> dgya(cmn.dgya, dimension(maxptn));
  arr_ref<double> dgza(cmn.dgza, dimension(maxptn));
  // COMMON ilist1
  arr_cref<int> next(cmn.next, dimension(maxptn));
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  // COMMON ilist5
  arr_ref<double> ct(cmn.ct, dimension(maxptn));
  arr_cref<double> ot(cmn.ot, dimension(maxptn));
  //
  // SAVE
  bool& allok = sve.allok;
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& icels0 = sve.icels0;
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  double& tm = sve.tm;
  double& tmin1 = sve.tmin1;
  //
  //C       this subroutine is used to update next(i), ct(i), ot(i),
  //C        and get tmin, nc for j
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist5/
  //C
  isco(cmn, i, j, allok, tm, t1, t2);
  //C
  if (allok) {
    //C       tm eq tmin, change nc to make sure fixtime get the collision with both
    //C       wall and particle
    //C
    if (tm < tmin) {
      tmin = tm;
      ct(j) = t2;
      nc = i;
      if (iconfg == 3 || iconfg == 5) {
        dgxa(j) = jxa * 10e0 * size1;
        dgya(j) = jya * 10e0 * size2;
        if (iconfg == 5) {
          dgza(j) = jza * 10e0 * size3;
        }
      }
    }
    //C
    if (tm <= ot(i)) {
      ct(i) = t1;
      icels0 = icels(i);
      i1 = icels0 / 10000;
      i2 = (icels0 - i1 * 10000) / 100;
      i3 = icels0 - i1 * 10000 - i2 * 100;
      wallc(cmn, i, i1, i2, i3, t, tmin1);
      fixtim(cmn, i, t, tmin1, tm, j);
      if (iconfg == 3 || iconfg == 5) {
        dgxa(i) = -jxa * 10e0 * size1;
        dgya(i) = -jya * 10e0 * size2;
        if (iconfg == 5) {
          dgza(i) = -jza * 10e0 * size3;
        }
      }
    }
    //C
    if (tm > ot(i) && next(i) == j) {
      ct(i) = t1;
      reor(cmn, t, tm, i, j);
    }
    //C
  }
  else if (next(i) == j) {
    //C
    tm = cmn.tlarge;
    //C
    reor(cmn, t, tm, i, j);
    //C
  }
  //C
}

struct chkcel_save
{
  int ick;
  int j;
  int jj;
  int jud2;
  int l;

  chkcel_save() :
    ick(fem::int0),
    j(fem::int0),
    jj(fem::int0),
    jud2(fem::int0),
    l(fem::int0)
  {}
};

void
chkcel(
  common& cmn,
  int const& il,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chkcel);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON ilist1
  const int maxptn = 400001;
  arr_cref<int> nic(cmn.nic, dimension(maxptn));
  // COMMON ilist2
  arr_cref<int, 3> icel(cmn.icel, dimension(10, 10, 10));
  //
  // SAVE
  int& ick = sve.ick;
  int& j = sve.j;
  int& jj = sve.jj;
  int& jud2 = sve.jud2;
  int& l = sve.l;
  //
  //C       this program is used to check through all the particles
  //C       in the cell (i1,i2,i3) and see if we can get a particle collision
  //C       with time less than the original input tmin ( the collision time of
  //C       il with the wall
  //C       and update the affected particles
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist4/
  //C
  if (iconfg == 3 || iconfg == 5) {
    jj = cmn.ichkpt;
    FEM_DO_SAFE(j, 1, jj) {
      ck(cmn, j, ick);
      //C     10/24/02 get rid of argument usage mismatch in ud2():
      jud2 = j;
      //C              if (ick .eq. 1) call ud2(j, il, t, tmin, nc)
      if (ick == 1) {
        ud2(cmn, jud2, il, t, tmin, nc);
      }
    }
    return;
  }
  //C
  if (i1 == 11 && i2 == 11 && i3 == 11) {
    l = cmn.icell;
  }
  else {
    l = icel(i1, i2, i3);
  }
  //C
  //C       if there is no particle
  if (l == 0) {
    return;
  }
  j = nic(l);
  //C       if there is only one particle
  if (j == 0) {
    ck(cmn, l, ick);
    if (ick == 1) {
      ud2(cmn, l, il, t, tmin, nc);
    }
    //C
    //C       if there are many particles
  }
  else {
    //C
    //C       we don't worry about the other colliding particle because it's
    //C       set in last(), and will be checked in ud2
    //C
    ck(cmn, l, ick);
    if (ick == 1) {
      ud2(cmn, l, il, t, tmin, nc);
    }
    //C
    while (j != l) {
      ck(cmn, j, ick);
      if (ick == 1) {
        ud2(cmn, j, il, t, tmin, nc);
      }
      j = nic(j);
    }
  }
  //C
}

struct chkout_save
{
  int i;
  int j;
  int k;
  int m1;
  int m2;
  int m3;

  chkout_save() :
    i(fem::int0),
    j(fem::int0),
    k(fem::int0),
    m1(fem::int0),
    m2(fem::int0),
    m3(fem::int0)
  {}
};

void
chkout(
  common& cmn,
  int const& l,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chkout);
  // SAVE
  int& i = sve.i;
  int& j = sve.j;
  int& k = sve.k;
  int& m1 = sve.m1;
  int& m2 = sve.m2;
  int& m3 = sve.m3;
  //
  //C       this subroutine is used to check the collisions with particles in
  //C       surface cells to see if we can get a smaller collision time than tmin
  //C       with particle nc, when the colliding particle is outside the cube
  //C       input l,t,tmin,nc
  //C       output tmin, nc
  //C
  //Cc      SAVE /prec2/
  //C
  m1 = 11;
  m2 = 11;
  m3 = 11;
  chkcel(cmn, l, m1, m2, m3, t, tmin, nc);
  //C
  FEM_DO_SAFE(i, 1, 10) {
    FEM_DO_SAFE(j, 1, 10) {
      FEM_DO_SAFE(k, 1, 10) {
        if (i == 1 || i == 10 || j == 1 || j == 10 || k == 1 || k == 10) {
          chkcel(cmn, l, i, j, k, t, tmin, nc);
        }
      }
    }
  }
  //C
}

struct chkin1_save
{
  int i;
  int itest;
  int j;
  int k;
  int m1;
  int m2;
  int m3;

  chkin1_save() :
    i(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0),
    m1(fem::int0),
    m2(fem::int0),
    m3(fem::int0)
  {}
};

void
chkin1(
  common& cmn,
  int const& l,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chkin1);
  // SAVE
  int& i = sve.i;
  int& itest = sve.itest;
  int& j = sve.j;
  int& k = sve.k;
  int& m1 = sve.m1;
  int& m2 = sve.m2;
  int& m3 = sve.m3;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube
  //C       and update the afftected particles through chkcel
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0;
  //C
  FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i >= 1 && i <= 10 && j >= 1 && j <= 10 && k >= 1 && k <= 10) {
          chkcel(cmn, l, i, j, k, t, tmin, nc);
        }
        else if (itest == 0) {
          m1 = 11;
          m2 = 11;
          m3 = 11;
          chkcel(cmn, l, m1, m2, m3, t, tmin, nc);
          itest = 1;
        }
      }
    }
  }
  //C
}

struct chkin2_save
{
  int i;
  int ia;
  int ib;
  int ic;
  int itest;
  int j;
  int k;

  chkin2_save() :
    i(fem::int0),
    ia(fem::int0),
    ib(fem::int0),
    ic(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
chkin2(
  common& cmn,
  int const& l,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chkin2);
  // SAVE
  int& i = sve.i;
  int& ia = sve.ia;
  int& ib = sve.ib;
  int& ic = sve.ic;
  int& j = sve.j;
  int& k = sve.k;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube
  //C       and update the afftected particles through chkcel
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  sve.itest = 0;
  //C
  FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        ia = i;
        ib = j;
        ic = k;
        if (k >= 1 && k <= 10) {
          if (i == 0) {
            ia = 10;
          }
          if (i == 11) {
            ia = 1;
          }
          if (j == 0) {
            ib = 10;
          }
          if (j == 11) {
            ib = 1;
          }
          chkcel(cmn, l, ia, ib, ic, t, tmin, nc);
        }
      }
    }
  }
  //C
}

struct chkin3_save
{
  int i;
  int ia;
  int ib;
  int ic;
  int itest;
  int j;
  int k;

  chkin3_save() :
    i(fem::int0),
    ia(fem::int0),
    ib(fem::int0),
    ic(fem::int0),
    itest(fem::int0),
    j(fem::int0),
    k(fem::int0)
  {}
};

void
chkin3(
  common& cmn,
  int const& l,
  int const& i1,
  int const& i2,
  int const& i3,
  double const& t,
  double& tmin,
  int& nc)
{
  FEM_CMN_SVE(chkin3);
  // SAVE
  int& i = sve.i;
  int& ia = sve.ia;
  int& ib = sve.ib;
  int& ic = sve.ic;
  int& j = sve.j;
  int& k = sve.k;
  //
  //C       this subroutine is used to check collisions for particle inside
  //C       the cube
  //C       and update the afftected particles through chkcel
  //C
  //C       itest is a flag to make sure the 111111 cell is checked only once
  sve.itest = 0;
  //C
  FEM_DO_SAFE(i, i1 - 1, i1 + 1) {
    FEM_DO_SAFE(j, i2 - 1, i2 + 1) {
      FEM_DO_SAFE(k, i3 - 1, i3 + 1) {
        if (i == 0) {
          ia = 10;
        }
        else if (i == 11) {
          ia = 1;
        }
        else {
          ia = i;
        }
        if (j == 0) {
          ib = 10;
        }
        else if (j == 11) {
          ib = 1;
        }
        else {
          ib = j;
        }
        if (k == 0) {
          ic = 10;
        }
        else if (k == 11) {
          ic = 1;
        }
        else {
          ic = k;
        }
        chkcel(cmn, l, ia, ib, ic, t, tmin, nc);
      }
    }
  }
  //C
}

struct ulist1_save
{
  int i1;
  int i2;
  int i3;
  int icels0;
  int k;
  int nc;
  double tmin;
  double tmin1;

  ulist1_save() :
    i1(fem::int0),
    i2(fem::int0),
    i3(fem::int0),
    icels0(fem::int0),
    k(fem::int0),
    nc(fem::int0),
    tmin(fem::double0),
    tmin1(fem::double0)
  {}
};

void
ulist1(
  common& cmn,
  int const& l,
  double const& t)
{
  FEM_CMN_SVE(ulist1);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON ilist1
  const int maxptn = 400001;
  arr_cref<int> icsta(cmn.icsta, dimension(maxptn));
  arr_cref<int> icels(cmn.icels, dimension(maxptn));
  //
  // SAVE
  int& i1 = sve.i1;
  int& i2 = sve.i2;
  int& i3 = sve.i3;
  int& icels0 = sve.icels0;
  int& nc = sve.nc;
  double& tmin = sve.tmin;
  double& tmin1 = sve.tmin1;
  //
  //C       this subroutine is used to update the interaction list when particle
  //C       l is disturbed.
  //C
  //Cc      SAVE /para5/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist5/
  //C
  icels0 = icels(l);
  i1 = icels0 / 10000;
  i2 = (icels0 - i1 * 10000) / 100;
  i3 = icels0 - i1 * 10000 - i2 * 100;
  //C       save collision info for use when the collision is a collision with wall
  //C       otherwise wallc will change icsta
  sve.k = fem::mod(icsta(l), 10);
  //C
  wallc(cmn, l, i1, i2, i3, t, tmin1);
  tmin = tmin1;
  nc = 0;
  //C
  if (i1 == 11 && i2 == 11 && i3 == 11) {
    chkout(cmn, l, t, tmin, nc);
  }
  else {
    if (iconfg == 1) {
      chkin1(cmn, l, i1, i2, i3, t, tmin, nc);
    }
    else if (iconfg == 2) {
      chkin2(cmn, l, i1, i2, i3, t, tmin, nc);
    }
    else if (iconfg == 4) {
      chkin3(cmn, l, i1, i2, i3, t, tmin, nc);
    }
    else if (iconfg == 3 || iconfg == 5) {
      chkcel(cmn, l, i1, i2, i3, t, tmin, nc);
    }
  }
  //C
  fixtim(cmn, l, t, tmin1, tmin, nc);
  //C
}

struct ulist_save
{
  int l;

  ulist_save() :
    l(fem::int0)
  {}
};

void
ulist(
  common& cmn,
  double const& t)
{
  FEM_CMN_SVE(ulist);
  // COMMON ilist1
  int& jscat = cmn.jscat;
  int& ictype = cmn.ictype;
  //
  // SAVE
  int& l = sve.l;
  //
  //C     this subroutine is used to update a new collision time list
  //C       notice this t has been updated
  //C
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist4/
  //C
  if (ictype == 1 || ictype == 2 || ictype == 5 || ictype == 6) {
    l = cmn.ifmpt;
    ulist1(cmn, l, t);
  }
  if (ictype != 1) {
    l = cmn.iscat;
    ulist1(cmn, l, t);
    if (jscat != 0) {
      l = jscat;
      ulist1(cmn, l, t);
    }
  }
  //C
}

struct zpcrun_save
{
  int iscat0;
  int jscat0;
  int niscat;
  int njscat;
  double t1;

  zpcrun_save() :
    iscat0(fem::int0),
    jscat0(fem::int0),
    niscat(fem::int0),
    njscat(fem::int0),
    t1(fem::double0)
  {}
};

//C
//C*****************************************************************************
//C
void
zpcrun(
  common& cmn,
  fem::star_type const& /* UNHANDLED: star argument */)
{
  FEM_CMN_SVE(zpcrun);
  common_write write(cmn);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON para7
  int& ioscar = cmn.ioscar;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> xmass(cmn.xmass, dimension(maxptn));
  arr_cref<int> ityp(cmn.ityp, dimension(maxptn));
  // COMMON ilist1
  int& iscat = cmn.iscat;
  int& jscat = cmn.jscat;
  arr_cref<int> next(cmn.next, dimension(maxptn));
  int& ictype = cmn.ictype;
  arr_cref<int> icsta(cmn.icsta, dimension(maxptn));
  // COMMON ilist4
  int& ifmpt = cmn.ifmpt;
  int& ichkpt = cmn.ichkpt;
  // COMMON ilist5
  arr_cref<double> ct(cmn.ct, dimension(maxptn));
  double& tlarge = cmn.tlarge;
  // COMMON ilist6
  double& t = cmn.t;
  int& iopern = cmn.iopern;
  int& icolln = cmn.icolln;
  //
  // SAVE
  int& iscat0 = sve.iscat0;
  int& jscat0 = sve.jscat0;
  int& niscat = sve.niscat;
  int& njscat = sve.njscat;
  double& t1 = sve.t1;
  //
  static const char* format_200 = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))";
  static const char* format_201 = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))";
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec4/
  //Cc      SAVE /prec5/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //Cc      SAVE /ilist6/
  //Cc      SAVE /ana1/
  //Cc      SAVE /anim/
  //C
  //C       save last collision info
  if (fem::mod(ictype, 2) == 0) {
    savrec(cmn, iscat);
    savrec(cmn, jscat);
  }
  //C
  //C1      get operation type
  getict(cmn, t1);
  //C2      check freezeout condition
  if (iconfg == 1 && t1 > tlarge / 2e0) {
    return;
  }
  const double tend1 = 250e0;
  if (iconfg == 2 || iconfg == 3) {
    if (t1 > tend1) {
      return;
    }
    //C           if (ichkpt .eq. mul) then
    //C              ii = 0
    //C              do i = 1, mul
    //C                 gztemp = gz(i) + vz(i) * (t1 - ft(i))
    //C                 if (sqrt(t1 ** 2 - gztemp ** 2) .lt. tend) then
    //C                    ii = 1
    //C                    goto 1000
    //C                 end if
    //C              end do
    //C 1000              continue
    //C              if (ii .eq. 0) return 1
    //C           end if
  }
  const double tend2 = 6.1e0;
  if (iconfg == 4 || iconfg == 5) {
    if (t1 > tend2) {
      return;
    }
  }
  //C
  //Clin-6/06/02 local freezeout for string melting,
  //C     decide what partons have frozen out at time t1:
  if (cmn.isoft == 5) {
    local(t1);
  }
  //C
  //C3      update iopern, t
  //C
  iopern++;
  t = t1;
  if (fem::mod(ictype, 2) == 0) {
    icolln++;
    //C
    //C     4/18/01-ctest off
    //C           write (2006, 1233) 'iscat=', iscat, 'jscat=', jscat,
    //C           write (2006, *) 'iscat=', iscat, ' jscat=', jscat,
    //C     1 ityp(iscat), ityp(jscat)
    //C           write (2006, 1233) 'iscat=', max(indx(iscat), indx(jscat)),
    //C     &        'jscat=', min(indx(iscat), indx(jscat))
    //C
    //C           write (2006, 1234) ' icolln=', icolln, 't=', t
    //C
    //C 1233           format (a10, i10, a10, i10)
    //C 1234           format (a15, i10, a5, f23.17, a5, f23.17)
  }
  //C
  //C4.1    deal with formation
  if (iconfg == 1 || iconfg == 2 || iconfg == 4) {
    if (ictype == 1 || ictype == 2 || ictype == 5 || ictype == 6) {
      celasn(cmn);
    }
  }
  //C
  //C4.2    deal with collisions
  //C
  if (ictype != 1) {
    //C
    iscat0 = iscat;
    jscat0 = jscat;
    //C
    //C        iscat is the larger one so that if it's a wall collision,
    //C       it's still ok
    iscat = fem::max0(iscat0, jscat0);
    jscat = fem::min0(iscat0, jscat0);
    //C
    //Ctest off check icsta(i): 0 with f77 compiler
    //C        write(9,*) 'BB:ictype,t1,iscat,jscat,icsta(i)=',
    //C     1 ictype,t1,iscat,jscat,icsta(iscat)
    //C
    //C       check collision time table error 'tterr'
    //Clin-4/2008 to avoid out-of-bound error in next():
    //C           if (jscat .ne. 0 .and. next(jscat) .ne. iscat)
    //C     &        then
    //C              print *, 'iscat=', iscat, 'jscat=', jscat,
    //C     &             'next(', jscat, ')=', next(jscat)
    //C
    //C              if (ct(iscat) .lt. tlarge / 2d0) stop 'tterr'
    //C              if (ct(jscat) .lt. tlarge / 2d0) stop 'tterr'
    //C           end if
    if (jscat != 0) {
      if (next(jscat) != iscat) {
        write(6, star), "iscat=", iscat, "jscat=", jscat, "next(",
          jscat, ")=", next(jscat);
        if (ct(iscat) < tlarge / 2e0) {
          FEM_STOP("tterr");
        }
        if (ct(jscat) < tlarge / 2e0) {
          FEM_STOP("tterr");
        }
      }
    }
    //Clin-4/2008-end
    //C
    //C4.2.1     collisions with wall
    //C
    //C     8/19/02 avoid actual argument in common blocks of cellre:
    niscat = iscat;
    njscat = jscat;
    //C           if (icsta(iscat) .ne. 0) call cellre(iscat, t)
    //C           if (jscat .ne. 0) then
    //C              if (icsta(jscat) .ne. 0) call cellre(jscat, t)
    //C           end if
    if (icsta(iscat) != 0) {
      cellre(cmn, niscat, t);
    }
    if (jscat != 0) {
      if (icsta(jscat) != 0) {
        cellre(cmn, njscat, t);
      }
    }
    //C
    //C4.2.2     collision between particles
    //C
    //Clin-6/2009 write out info for each collision:
    //C           if (mod(ictype, 2) .eq. 0) call scat(t, iscat, jscat)
    if (fem::mod(ictype, 2) == 0) {
      if (ioscar == 3) {
        write(95, star), "event,miss,iscat,jscat=", cmn.iaevt,
          cmn.miss, iscat, jscat;
        if (fem::dmax1(fem::abs(gx(iscat)), fem::abs(gy(iscat)),
            fem::abs(gz(iscat)), fem::abs(ft(iscat)), fem::abs(gx(jscat)),
            fem::abs(gy(jscat)), fem::abs(gz(jscat)), fem::abs(ft(
            jscat))) < 9999) {
          write(95, format_200), ityp(iscat), px(iscat), py(iscat),
            pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat),
            ft(iscat);
          write(95, format_200), ityp(jscat), px(jscat), py(jscat),
            pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat),
            ft(jscat);
        }
        else {
          write(95, format_201), ityp(iscat), px(iscat), py(iscat),
            pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat),
            ft(iscat);
          write(95, format_201), ityp(jscat), px(jscat), py(jscat),
            pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat),
            ft(jscat);
        }
      }
      //C
      scat(cmn, t, iscat, jscat);
      //C
      if (ioscar == 3) {
        if (fem::dmax1(fem::abs(gx(iscat)), fem::abs(gy(iscat)),
            fem::abs(gz(iscat)), fem::abs(ft(iscat)), fem::abs(gx(jscat)),
            fem::abs(gy(jscat)), fem::abs(gz(jscat)), fem::abs(ft(
            jscat))) < 9999) {
          write(95, format_200), ityp(iscat), px(iscat), py(iscat),
            pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat),
            ft(iscat);
          write(95, format_200), ityp(jscat), px(jscat), py(jscat),
            pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat),
            ft(jscat);
        }
        else {
          write(95, format_201), ityp(iscat), px(iscat), py(iscat),
            pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat),
            ft(iscat);
          write(95, format_201), ityp(jscat), px(jscat), py(jscat),
            pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat),
            ft(jscat);
        }
      }
    }
    //C
  }
  //C
  //C5      update the interaction list
  ulist(cmn, t);
  //C
  //C6      update ifmpt. ichkpt
  //C       old ichkpt and ifmpt are more conveniently used in ulist
  if (ifmpt <= cmn.mul) {
    if (ictype != 0 && ictype != 3 && ictype != 4) {
      ichkpt++;
      ifmpt++;
    }
  }
  //C
}

void
zpca1b(
  common& cmn,
  double const& rapi,
  double const& et,
  int const& ian)
{
  // COMMON para6
  double& centy = cmn.centy;
  // COMMON ana2
  arr_ref<double> det(cmn.det, dimension(12));
  arr_ref<double> dn(cmn.dn, dimension(12));
  arr_ref<double> det1(cmn.det1, dimension(12));
  arr_ref<double> dn1(cmn.dn1, dimension(12));
  arr_ref<double> det2(cmn.det2, dimension(12));
  arr_ref<double> dn2(cmn.dn2, dimension(12));
  //
  //C
  //Cc      SAVE /para6/
  //Cc      SAVE /ilist6/
  //Cc      SAVE /ana2/
  //C
  if (rapi > centy - 0.5e0 && rapi < centy + 0.5e0) {
    det2(ian) += et;
    dn2(ian) += 1e0;
    //Cdtrans
    if (ian == 10) {
      //Cd              write (10, *) t, det2(ian)
    }
    if (ian == 11) {
      //Cd              write (11, *) t, det2(ian)
    }
    if (ian == 12) {
      //Cd              write (12, *) t, det2(ian)
    }
    //Cdtransend
    if (rapi > centy - 0.25e0 && rapi < centy + 0.25e0) {
      det1(ian) += et;
      dn1(ian) += 1e0;
      if (rapi > centy - 0.1e0 && rapi < centy + 0.1e0) {
        det(ian) += et;
        dn(ian) += 1e0;
      }
    }
  }
  //C
}

struct zpca1c_save
{
  arr<double> en;
  int i;
  int j;

  zpca1c_save() :
    en(dimension(4), fem::fill0),
    i(fem::int0),
    j(fem::int0)
  {}
};

void
zpca1c(
  common& cmn,
  double const& p0,
  double const& p1,
  double const& p2,
  double const& p3,
  int const& ian)
{
  FEM_CMN_SVE(zpca1c);
  // COMMON ana3
  arr_ref<double, 3> em(cmn.em, dimension(4, 4, 12));
  //
  // SAVE
  arr_ref<double> en(sve.en, dimension(4));
  int& i = sve.i;
  int& j = sve.j;
  //
  //C
  //Cc      SAVE /ana3/
  //C
  en(1) = p0;
  en(2) = p1;
  en(3) = p2;
  en(4) = p3;
  //C
  FEM_DO_SAFE(i, 1, 4) {
    FEM_DO_SAFE(j, 1, 4) {
      em(i, j, ian) += en(i) * en(j) / p0;
    }
  }
  //C
}

struct zpca1a_save
{
  double et;
  int ian;
  int ipic;
  double p0;
  double p1;
  double p2;
  double p3;
  double rapi;
  double t1;
  double t2;

  zpca1a_save() :
    et(fem::double0),
    ian(fem::int0),
    ipic(fem::int0),
    p0(fem::double0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    rapi(fem::double0),
    t1(fem::double0),
    t2(fem::double0)
  {}
};

void
zpca1a(
  common& cmn,
  int const& i)
{
  FEM_CMN_SVE(zpca1a);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  // COMMON prec3
  arr_cref<double> fts(cmn.fts, dimension(maxptn));
  arr_cref<double> pxs(cmn.pxs, dimension(maxptn));
  arr_cref<double> pys(cmn.pys, dimension(maxptn));
  arr_cref<double> pzs(cmn.pzs, dimension(maxptn));
  arr_cref<double> es(cmn.es, dimension(maxptn));
  // COMMON prec5
  arr_cref<double> tau(cmn.tau, dimension(maxptn));
  // COMMON prec6
  arr_cref<double> raps(cmn.raps, dimension(maxptn));
  arr_cref<double> taus(cmn.taus, dimension(maxptn));
  // COMMON ana1
  arr_cref<double> ts(cmn.ts, dimension(12));
  //
  // SAVE
  double& et = sve.et;
  int& ian = sve.ian;
  int& ipic = sve.ipic;
  double& p0 = sve.p0;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& rapi = sve.rapi;
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  //
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para5/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec3/
  //Cc      SAVE /prec5/
  //Cc      SAVE /prec6/
  //Cc      SAVE /ana1/
  //C
  if (iconfg == 1) {
    t1 = fts(i);
    t2 = ft(i);
    ipic = 11;
  }
  else if (iconfg == 2 || iconfg == 3) {
    //Cd           t1 = fts(i)
    //Cd           t2 = ft(i)
    t1 = taus(i);
    t2 = tau(i);
    ipic = 12;
  }
  else if (iconfg == 4 || iconfg == 5) {
    t1 = fts(i);
    t2 = ft(i);
    ipic = 12;
  }
  //C
  if (iconfg <= 3) {
    FEM_DO_SAFE(ian, 1, ipic) {
      if (t1 <= ts(ian) && t2 > ts(ian)) {
        rapi = raps(i);
        //C     7/20/01:
        //C                 et = sqrt(pxs(i) ** 2 + pys(i) ** 2 + xmp ** 2)
        et = fem::dsqrt(fem::pow2(pxs(i)) + fem::pow2(pys(i)) +
          fem::pow2(cmn.xmp));
        zpca1b(cmn, rapi, et, ian);
      }
    }
  }
  else {
    FEM_DO_SAFE(ian, 1, ipic) {
      if (t1 <= ts(ian) && t2 > ts(ian)) {
        p0 = es(i);
        p1 = pxs(i);
        p2 = pys(i);
        p3 = pzs(i);
        zpca1c(cmn, p0, p1, p2, p3, ian);
      }
    }
  }
  //C
}

//C
//C*****************************************************************************
//C
void
zpca1(
  common& cmn)
{
  //C
  //Cc      SAVE /ilist1/
  //C
  if (fem::mod(cmn.ictype, 2) == 0) {
    zpca1a(cmn, cmn.iscat);
    zpca1a(cmn, cmn.jscat);
    //Clin-5/2009 ctest off v2 for parton:
    //C           call flowp(1)
  }
  //C
}

struct zpca2a_save
{
  double et;
  int i;
  int ian;
  int ipic;
  int j;
  double rapi;
  double t1;
  double t2;

  zpca2a_save() :
    et(fem::double0),
    i(fem::int0),
    ian(fem::int0),
    ipic(fem::int0),
    j(fem::int0),
    rapi(fem::double0),
    t1(fem::double0),
    t2(fem::double0)
  {}
};

void
zpca2a(
  common& cmn)
{
  FEM_CMN_SVE(zpca2a);
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON para6
  double& centy = cmn.centy;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  // COMMON prec5
  arr_cref<double> rap(cmn.rap, dimension(maxptn));
  arr_cref<double> tau(cmn.tau, dimension(maxptn));
  // COMMON ilist5
  double& tlarge = cmn.tlarge;
  // COMMON ana1
  arr_cref<double> ts(cmn.ts, dimension(12));
  // COMMON ana2
  arr_cref<double> det(cmn.det, dimension(12));
  arr_cref<double> dn(cmn.dn, dimension(12));
  arr_ref<double> detdy(cmn.detdy, dimension(12));
  arr_ref<double> detdn(cmn.detdn, dimension(12));
  arr_ref<double> dndy(cmn.dndy, dimension(12));
  arr_cref<double> det1(cmn.det1, dimension(12));
  arr_cref<double> dn1(cmn.dn1, dimension(12));
  arr_ref<double> detdy1(cmn.detdy1, dimension(12));
  arr_ref<double> detdn1(cmn.detdn1, dimension(12));
  arr_ref<double> dndy1(cmn.dndy1, dimension(12));
  arr_cref<double> det2(cmn.det2, dimension(12));
  arr_cref<double> dn2(cmn.dn2, dimension(12));
  arr_ref<double> detdy2(cmn.detdy2, dimension(12));
  arr_ref<double> detdn2(cmn.detdn2, dimension(12));
  arr_ref<double> dndy2(cmn.dndy2, dimension(12));
  // COMMON ana4
  arr_ref<double> fdetdy(cmn.fdetdy, dimension(24));
  arr_ref<double> fdndy(cmn.fdndy, dimension(24));
  arr_ref<double> fdndpt(cmn.fdndpt, dimension(12));
  //
  // SAVE
  double& et = sve.et;
  int& i = sve.i;
  int& ian = sve.ian;
  int& ipic = sve.ipic;
  int& j = sve.j;
  double& rapi = sve.rapi;
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  //
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para2/
  //Cc      SAVE /para3/
  //Cc      SAVE /para5/
  //Cc      SAVE /para6/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec5/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //Cc      SAVE /ilist6/
  //Cc      SAVE /rndm1/
  //Cc      SAVE /rndm2/
  //Cc      SAVE /rndm3/
  //Cc      SAVE /ana1/
  //Cc      SAVE /ana2/
  //Cc      SAVE /ana4/
  //C
  FEM_DO_SAFE(i, 1, cmn.ichkpt) {
    rapi = rap(i);
    //C     7/20/01:
    //C           et = sqrt(px(i) ** 2 + py(i) ** 2 + xmp ** 2)
    et = fem::dsqrt(fem::pow2(px(i)) + fem::pow2(py(i)) + fem::pow2(cmn.xmp));
    //C
    FEM_DO_SAFE(j, 1, 24) {
      if (rapi > j + centy - 13e0 && rapi < j + centy - 12e0) {
        fdetdy(j) += et;
        fdndy(j) += 1e0;
      }
    }
    //C
    FEM_DO_SAFE(j, 1, 12) {
      if (et > 0.5e0 * (j - 1) && et < 0.5e0 * j) {
        fdndpt(j) += 1e0;
      }
    }
    //C
    if (iconfg == 1) {
      t1 = ft(i);
      t2 = tlarge;
      ipic = 11;
    }
    else {
      t1 = tau(i);
      t2 = tlarge;
      ipic = 12;
    }
    //C
    FEM_DO_SAFE(ian, 1, ipic) {
      if (t1 <= ts(ian) && t2 > ts(ian)) {
        zpca1b(cmn, rapi, et, ian);
      }
    }
    //C
    if (iconfg == 1) {
      zpca1b(cmn, rapi, et, 12);
    }
  }
  //C
  FEM_DO_SAFE(ian, 1, 12) {
    if (dn(ian) == 0e0 || dn1(ian) == 0e0 || dn2(ian) == 0e0) {
      //Clin-9/2012 suppress output:
      //C              print *, 'event=', ievt
      //C              print *, 'dn(', ian, ')=', dn(ian), 'dn1(', ian,
      //C     &           ')=', dn1(ian), 'dn2(', ian, ')=', dn2(ian)
    }
    detdy(ian) += det(ian);
    if (dn(ian) != 0) {
      detdn(ian) += det(ian) / dn(ian);
    }
    dndy(ian) += dn(ian);
    detdy1(ian) += det1(ian);
    if (dn1(ian) != 0) {
      detdn1(ian) += det1(ian) / dn1(ian);
    }
    dndy1(ian) += dn1(ian);
    detdy2(ian) += det2(ian);
    if (dn2(ian) != 0) {
      detdn2(ian) += det2(ian) / dn2(ian);
    }
    dndy2(ian) += dn2(ian);
  }
  //C
}

struct zpca2b_save
{
  int i;
  int ian;
  int ipic;
  double p0;
  double p1;
  double p2;
  double p3;
  double t1;
  double t2;

  zpca2b_save() :
    i(fem::int0),
    ian(fem::int0),
    ipic(fem::int0),
    p0(fem::double0),
    p1(fem::double0),
    p2(fem::double0),
    p3(fem::double0),
    t1(fem::double0),
    t2(fem::double0)
  {}
};

void
zpca2b(
  common& cmn)
{
  FEM_CMN_SVE(zpca2b);
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  // COMMON ana1
  arr_cref<double> ts(cmn.ts, dimension(12));
  //
  // SAVE
  int& i = sve.i;
  int& ian = sve.ian;
  int& ipic = sve.ipic;
  double& p0 = sve.p0;
  double& p1 = sve.p1;
  double& p2 = sve.p2;
  double& p3 = sve.p3;
  double& t1 = sve.t1;
  double& t2 = sve.t2;
  //
  //C
  //Cc      SAVE /prec2/
  //Cc      SAVE /ilist4/
  //Cc      SAVE /ilist5/
  //Cc      SAVE /ana1/
  //C
  FEM_DO_SAFE(i, 1, cmn.ichkpt) {
    t1 = ft(i);
    t2 = cmn.tlarge;
    ipic = 12;
    //C
    FEM_DO_SAFE(ian, 1, ipic) {
      if (t1 <= ts(ian) && t2 > ts(ian)) {
        p0 = e(i);
        p1 = px(i);
        p2 = py(i);
        p3 = pz(i);
        zpca1c(cmn, p0, p1, p2, p3, ian);
      }
    }
  }
  //C
}

struct zpca2c_save
{
  int aproj;
  int atarg;
  double bimp;
  fem::str<8> code;
  double ebeam;
  int event;
  int i;
  int nff;
  int ntestp;
  double phi;
  fem::str<4> reffra;
  fem::str<8> versn;
  int zproj;
  int ztarg;

  zpca2c_save() :
    aproj(fem::int0),
    atarg(fem::int0),
    bimp(fem::double0),
    code(fem::char0),
    ebeam(fem::double0),
    event(fem::int0),
    i(fem::int0),
    nff(fem::int0),
    ntestp(fem::int0),
    phi(fem::double0),
    reffra(fem::char0),
    versn(fem::char0),
    zproj(fem::int0),
    ztarg(fem::int0)
  {}
};

void
zpca2c(
  common& cmn)
{
  FEM_CMN_SVE(zpca2c);
  common_write write(cmn);
  // COMMON para1
  int& mul = cmn.mul;
  // COMMON prec2
  const int maxptn = 400001;
  arr_cref<double> gx(cmn.gx, dimension(maxptn));
  arr_cref<double> gy(cmn.gy, dimension(maxptn));
  arr_cref<double> gz(cmn.gz, dimension(maxptn));
  arr_cref<double> ft(cmn.ft, dimension(maxptn));
  arr_cref<double> px(cmn.px, dimension(maxptn));
  arr_cref<double> py(cmn.py, dimension(maxptn));
  arr_cref<double> pz(cmn.pz, dimension(maxptn));
  arr_cref<double> e(cmn.e, dimension(maxptn));
  arr_cref<double> xmass(cmn.xmass, dimension(maxptn));
  arr_cref<int> ityp(cmn.ityp, dimension(maxptn));
  //
  // SAVE
  int& aproj = sve.aproj;
  int& atarg = sve.atarg;
  double& bimp = sve.bimp;
  fem::str<8>& code = sve.code;
  double& ebeam = sve.ebeam;
  int& event = sve.event;
  int& i = sve.i;
  int& nff = sve.nff;
  int& ntestp = sve.ntestp;
  double& phi = sve.phi;
  fem::str<4>& reffra = sve.reffra;
  fem::str<8>& versn = sve.versn;
  int& zproj = sve.zproj;
  int& ztarg = sve.ztarg;
  //
  if (is_called_first_time) {
    nff = 0;
  }
  static const char* format_101 = "(a12)";
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /prec2/
  //C
  //C       file header
  if (nff == 0) {
    write(26, format_101), "OSCAR1997A";
    write(26, format_101), "final_id_p_x";
    code = "ZPC";
    versn = "1.0.1";
    aproj = -1;
    zproj = -1;
    atarg = -1;
    ztarg = -1;
    reffra = "cm";
    ebeam = 0e0;
    ntestp = 1;
    write(26,
      "(2(a8,2x),'(',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,e10.4,2x,i8)"),
      code, versn, aproj, zproj, atarg, ztarg, reffra, ebeam, ntestp;
    nff = 1;
    event = 1;
    bimp = 0e0;
    phi = 0e0;
  }
  //C
  //C       comment
  //C
  //C       event header
  write(26, "(i10,2x,i10,2x,f8.3,2x,f8.3)"), event, mul, bimp, phi;
  //C
  //C       particles
  FEM_DO_SAFE(i, 1, mul) {
    write(26, "(i10,2x,i10,2x,9(e12.6,2x))"), i, ityp(i), px(i), py(i),
      pz(i), e(i), xmass(i), gx(i), gy(i), gz(i), ft(i);
  }
  //C
  event++;
  //C
}

//C
//C*****************************************************************************
//C
void
zpca2(
  common& cmn)
{
  common_write write(cmn);
  //C
  //Cc      SAVE /para3/
  //Cc      SAVE /para5/
  //Cc      SAVE /para7/
  //Cc      SAVE /ilist6/
  //Cc      SAVE /rndm1/
  //Cc      SAVE /rndm2/
  //Cc      SAVE /rndm3/
  //Cc      SAVE /AREVT/
  //C
  if (cmn.iconfg <= 3) {
    zpca2a(cmn);
  }
  else {
    zpca2b(cmn);
  }
  //C
  if (cmn.ioscar == 1) {
    zpca2c(cmn);
  }
  //C
  //Cbzdbg2/17/99
  //C        write (25, *) 'Event', nsevt - 1 + ievt,
  //C    &         ', run', isbrun,
  //C        WRITE (25, *) ' Event ', IAEVT, ', run ', IARUN,
  //C     &     ',\n\t number of operations = ', iopern,
  //C     &     ',\n\t number of collisions between particles = ',
  //C     &         icolln,
  //C     &     ',\n\t freezeout time=', t,
  //C     &     ',\n\t ending at the ', number, 'th random number',
  //C     &     ',\n\t ending collision iff=', iff
  write(25, star), " Event ", cmn.iaevt, ", run ", cmn.iarun;
  write(25, star), "    number of operations = ", cmn.iopern;
  write(25, star), "    number of collisions between particles = ", cmn.icolln;
  write(25, star), "    freezeout time=", cmn.t;
  write(25, star), "    ending at the ", cmn.number, "th random number";
  write(25, star), "    ending collision iff=", cmn.iff;
  //C
}

struct zpcou1_save
{
  double dpt;
  double dy;
  double dy1;
  double dy2;
  int ntotal;

  zpcou1_save() :
    dpt(fem::double0),
    dy(fem::double0),
    dy1(fem::double0),
    dy2(fem::double0),
    ntotal(fem::int0)
  {}
};

void
zpcou1(
  common& cmn)
{
  FEM_CMN_SVE(zpcou1);
  //C
  //Cc      SAVE /para3/
  //Cc      SAVE /ana1/
  //Cc      SAVE /ana2/
  //Cc      SAVE /ana4/
  //C
  sve.dpt = 0.5e0;
  sve.dy2 = 1e0;
  sve.dy1 = 0.5e0;
  sve.dy = 0.2e0;
  sve.ntotal = cmn.nevnt * cmn.nsbrun;
  //C
}

struct zpcou2_save
{
  int i;
  int ian;
  int ntotal;
  double vol;

  zpcou2_save() :
    i(fem::int0),
    ian(fem::int0),
    ntotal(fem::int0),
    vol(fem::double0)
  {}
};

void
zpcou2(
  common& cmn)
{
  FEM_CMN_SVE(zpcou2);
  common_write write(cmn);
  // COMMON ana1
  arr_cref<double> ts(cmn.ts, dimension(12));
  // COMMON ana3
  arr_cref<double, 3> em(cmn.em, dimension(4, 4, 12));
  //
  // SAVE
  int& i = sve.i;
  int& ian = sve.ian;
  int& ntotal = sve.ntotal;
  double& vol = sve.vol;
  //
  //C
  //Cc      SAVE /para3/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ana1/
  //Cc      SAVE /ana3/
  //C
  cmn.io.open(28, "ana4/em.dat")
    .status("unknown");
  vol = 1000.e0 * cmn.size1 * cmn.size2 * cmn.size3;
  ntotal = cmn.nevnt * cmn.nsbrun;
  //C
  FEM_DO_SAFE(ian, 1, 12) {
    write(28, star), "*** for time ", ts(ian), "fm(s)";
    FEM_DO_SAFE(i, 1, 4) {
      write(28, star), em(i, 1, ian) / vol / ntotal, em(i, 2, ian) /
        vol / ntotal, em(i, 3, ian) / vol / ntotal, em(i, 4, ian) /
        vol / ntotal;
    }
  }
  //C
}

//C
//C*****************************************************************************
//C
void
zpcou(
  common& cmn)
{
  //C
  //Cc      SAVE /para5/
  //C
  if (cmn.iconfg <= 3) {
    zpcou1(cmn);
  }
  else {
    zpcou2(cmn);
  }
  //C
}

struct zpcmn_save
{
  int i;
  int j;

  zpcmn_save() :
    i(fem::int0),
    j(fem::int0)
  {}
};

//C.................... zpc.f
//C        PROGRAM ZPC
void
zpcmn(
  common& cmn)
{
  FEM_CMN_SVE(zpcmn);
  int& i = sve.i;
  int& j = sve.j;
  //C       Version: 1.0.1
  //C       Author: Bin Zhang
  //C       (suggestions, problems -> bzhang@nt1.phys.columbia.edu)
  //Clin-4/20/01        PARAMETER (NMAXGL = 16000)
  //Cc      SAVE /para3/
  //C
  //C       loop over events
  FEM_DO_SAFE(i, 1, cmn.nevnt) {
    cmn.ievt = i;
    //C       generation of the initial condition for one event
    inievt(cmn);
    //C      loop over many runs of the same event
    FEM_DO_SAFE(j, 1, cmn.nsbrun) {
      cmn.isbrun = j;
      //C       initialization for one run of an event
      inirun(cmn);
      //Clin-4/2008 not used:
      //C             CALL HJAN1A
      statement_3000:
      //C       do one collision
      zpcrun(cmn, star /* 4000 UNHANDLED */);
      zpca1(cmn);
      goto statement_3000;
      zpca2(cmn);
    }
  }
  zpcou(cmn);
  //Clin-5/2009 ctest off
  //C     5/17/01 calculate v2 for parton already frozen out:
  //C        call flowp(3)
  //C.....to get average values for different strings
  zpstrg();
}

struct blockdata_zpcbdt_save
{
};

//C
//C*****************************************************************************
//C
void
blockdata_zpcbdt(
  common& cmn)
{
  FEM_CMN_SVE(blockdata_zpcbdt);
  // COMMON ana1
  arr_ref<double> ts(cmn.ts, dimension(12));
  //
  if (is_called_first_time) {
    cmn.centy = 0e0;
    cmn.number = 0;
    {
      static const double values[] = {
        0.11e0, 0.12e0, 0.15e0, 0.2e0, 0.3e0, 0.4e0, 0.6e0, 0.8e0,
          1e0, 2e0, 4e0, 6e0
      };
      fem::data_of_type<double>(FEM_VALUES_AND_SIZE),
        ts;
    }
  }
  //C       set initial values in block data
  //C
  //Cc      SAVE /para1/
  //Cc      SAVE /para2/
  //Cc      SAVE /para3/
  //Cc      SAVE /para4/
  //Cc      SAVE /para5/
  //Cc      SAVE /para6/
  //Clin-6/2009 nsmbbbar and nsmmeson respectively give the total number of
  //C     baryons/anti-baryons and mesons for each event:
  //C        common /para7/ ioscar
  //Cc      SAVE /para7/
  //Cc      SAVE /prec1/
  //Cc      SAVE /prec2/
  //Cc      SAVE /prec3/
  //Cc      SAVE /prec4/
  //Cc      SAVE /prec5/
  //Cc      SAVE /prec6/
  //Cc      SAVE /aurec1/
  //Cc      SAVE /aurec2/
  //Cc      SAVE /ilist1/
  //Cc      SAVE /ilist2/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /ilist4/
  //C     6/07/02 initialize in ftime to expedite compiling:
  //C        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
  //Cc      SAVE /ilist5/
  //Cc      SAVE /ilist6/
  //Cc      SAVE /ilist7/
  //Cc      SAVE /ilist8/
  //Cc      SAVE /rndm1/
  //Cc      SAVE /rndm2/
  //Cc      SAVE /rndm3/
  //Cc      SAVE /ana1/
  //Cc      SAVE /ana2/
  //Cc      SAVE /ana3/
  //Cc      SAVE /ana4/
  //C     6/07/02 initialize in ftime to expedite compiling:
  //C        data (ct(i), i = 1, MAXPTN)/MAXPTN*0d0/
  //C        data (ot(i), i = 1, MAXPTN)/MAXPTN*0d0/
  //C        data tlarge/1000000.d0/
  //C
}

struct readpa_save
{
  double a;
  int i;
  int irused;
  int isedng;
  int iseed;
  int iseed2;
  fem::str<50> str;

  readpa_save() :
    a(fem::double0),
    i(fem::int0),
    irused(fem::int0),
    isedng(fem::int0),
    iseed(fem::int0),
    iseed2(fem::int0),
    str(fem::char0)
  {}
};

void
readpa(
  common& cmn)
{
  FEM_CMN_SVE(readpa);
  common_write write(cmn);
  // COMMON para4
  int& ireflg = cmn.ireflg;
  // COMMON para5
  int& iconfg = cmn.iconfg;
  // COMMON ilist3
  double& size1 = cmn.size1;
  double& size2 = cmn.size2;
  double& size3 = cmn.size3;
  double& v1 = cmn.v1;
  double& v2 = cmn.v2;
  double& v3 = cmn.v3;
  //
  // SAVE
  double& a = sve.a;
  int& irused = sve.irused;
  int& isedng = sve.isedng;
  int& iseed = sve.iseed;
  int& iseed2 = sve.iseed2;
  //
  //C
  //Cc      SAVE /para2/
  //Cc      SAVE /para3/
  //Cc      SAVE /para4/
  //Cc      SAVE /para5/
  //Cc      SAVE /para7/
  //Cc      SAVE /ilist3/
  //Cc      SAVE /rndm1/
  //Cc      SAVE /rndm2/
  //Cc      SAVE /rndm3/
  //C
  iseed = cmn.iseedp;
  //C       this is the initialization file containing the initial values of
  //C          the parameters
  //Cbz1/31/99
  //C        open (5, file = 'zpc.ini', status = 'unknown')
  //Cbz1/31/99end
  //C
  //C       this is the final data file containing general info about the cascade
  //Cbz1/31/99
  //C        open (6, file = 'zpc.res', status = 'unknown')
  cmn.io.open(25, "ana/zpc.res")
    .status("unknown");
  //Cbz1/31/99end
  //C
  //C       this is the input file containing initial particle records
  //Cbz1/25/99
  //C        open (7, file = 'zpc.inp', status = 'unknown')
  //Cbz1/25/99end
  //C
  //C       this gives the optional OSCAR standard output
  //Cbz1/31/99
  //C        open (8, file = 'zpc.oscar', status = 'unknown')
  if (cmn.ioscar == 1) {
    cmn.io.open(26, "ana/parton.oscar")
      .status("unknown");
    cmn.io.open(19, "ana/hadron.oscar")
      .status("unknown");
  }
  //Cbz1/31/99end
  //C
  //C     2/11/03 combine zpc initialization into ampt.ini:
  //C        open (29, file = 'zpc.ini', status = 'unknown')
  //C        read (29, *) str, xmp
  cmn.xmp = 0e0;
  //C        read (29, *) str, xmu
  //C        read (29, *) str, alpha
  cmn.cutof2 = 4.5e0 * fem::pow2((cmn.alpha / cmn.xmu));
  //C        read (29, *) str, rscut2
  cmn.rscut2 = 0.01e0;
  //C        read (29, *) str, nsevt
  cmn.nsevt = 1;
  //C        read (29, *) str, nevnt
  cmn.nevnt = 1;
  //C        read (29, *) str, nsbrun
  cmn.nsbrun = 1;
  //C        read (29, *) str, iftflg
  cmn.iftflg = 0;
  //C        read (29, *) str, ireflg
  ireflg = 1;
  //Cbz1/31/99
  if (ireflg == 0) {
    cmn.io.open(27, "zpc.inp")
      .status("UNKNOWN");
  }
  //Cbz1/31/99end
  //C        read (29, *) str, igeflg
  cmn.igeflg = 0;
  //C        read (29, *) str, ibstfg
  cmn.ibstfg = 0;
  //C        read (29, *) str, iconfg
  iconfg = 1;
  //C        read (29, *) str, iordsc
  cmn.iordsc = 11;
  //C        read (29, *) str, ioscar
  //C        read (29, *) str, v1, v2, v3
  v1 = 0.2e0;
  v2 = 0.2e0;
  v3 = 0.2e0;
  //C        read (29, *) str, size1, size2, size3
  size1 = 1.5e0;
  size2 = 1.5e0;
  size3 = 0.7e0;
  if (size1 == 0e0 || size2 == 0e0 || size3 == 0e0) {
    if (size1 != 0e0 || size2 != 0e0 || size3 != 0e0 || v1 != 0e0 ||
        v2 != 0e0 || v3 != 0e0) {
      write(6, star), "to get rid of space division:";
      write(6, star), "set all sizes and vs to 0";
      FEM_STOP("chker");
    }
  }
  cmn.size = fem::min(size1, size2, size3);
  //C        read (29, *) str, iff
  cmn.iff = -1;
  //C        read (29, *) str, iseed
  //C
  //C     10/24/02 get rid of argument usage mismatch in ran1():
  isedng = -iseed;
  //C        a = ran1(-iseed)
  a = ran1(cmn, isedng);
  //C        read (29, *) str, irused
  irused = 2;
  FEM_DO_SAFE(sve.i, 1, irused - 1) {
    //C           a = ran1(2)
    iseed2 = 2;
    a = ran1(cmn, iseed2);
  }
  //C     10/24/02-end
  //C
  if (iconfg == 2 || iconfg == 3) {
    v1 = 0e0;
    v2 = 0e0;
  }
  //C
  if (iconfg == 4 || iconfg == 5) {
    v1 = 0e0;
    v2 = 0e0;
    v3 = 0e0;
  }
  //C
  cmn.io.close(5);
  //C
}

void
inipar(
  common& cmn)
{
  //C
  //Cc      SAVE /para4/
  //Cc      SAVE /para6/
  //C
  if (cmn.ibstfg != 0) {
    cmn.centy = -6e0;
  }
  //C
}

struct inian1_save
{
  double a;
  int i;

  inian1_save() :
    a(fem::double0),
    i(fem::int0)
  {}
};

void
inian1(
  common& cmn)
{
  FEM_CMN_SVE(inian1);
  // COMMON ana1
  arr_ref<double> ts(cmn.ts, dimension(12));
  //
  // SAVE
  double& a = sve.a;
  int& i = sve.i;
  //
  //C
  //Cc      SAVE /para4/
  //Cc      SAVE /ana1/
  if (cmn.ibstfg != 0) {
    a = fem::cosh(6e0);
    FEM_DO_SAFE(i, 1, 12) {
      ts(i) = ts(i) * a;
    }
  }
  //C
}

//C
//C*****************************************************************************
//C
void
inizpc(
  common& cmn)
{
  //C
  readpa(cmn);
  //C
  inipar(cmn);
  //C
  inian1(cmn);
  //C
}

} // namespace AMPT
