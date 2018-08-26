#include <fem.hpp> // Fortran EMulation library of fable module

namespace AMPT {

using namespace fem::major_types;

void
ardata(...)
{
  throw std::runtime_error(
    "Missing function implementation: ardata");
}

void
arini(...)
{
  throw std::runtime_error(
    "Missing function implementation: arini");
}

void
arini2(...)
{
  throw std::runtime_error(
    "Missing function implementation: arini2");
}

void
artan1(...)
{
  throw std::runtime_error(
    "Missing function implementation: artan1");
}

void
artan2(...)
{
  throw std::runtime_error(
    "Missing function implementation: artan2");
}

void
artmn(...)
{
  throw std::runtime_error(
    "Missing function implementation: artmn");
}

void
artout(...)
{
  throw std::runtime_error(
    "Missing function implementation: artout");
}

void
artset(...)
{
  throw std::runtime_error(
    "Missing function implementation: artset");
}

void
getnp(...)
{
  throw std::runtime_error(
    "Missing function implementation: getnp");
}

void
hidata(...)
{
  throw std::runtime_error(
    "Missing function implementation: hidata");
}

void
hijing(...)
{
  throw std::runtime_error(
    "Missing function implementation: hijing");
}

void
hijset(...)
{
  throw std::runtime_error(
    "Missing function implementation: hijset");
}

void
inizpc(...)
{
  throw std::runtime_error(
    "Missing function implementation: inizpc");
}

void
ludata(...)
{
  throw std::runtime_error(
    "Missing function implementation: ludata");
}

void
ppbdat(...)
{
  throw std::runtime_error(
    "Missing function implementation: ppbdat");
}

void
pydata(...)
{
  throw std::runtime_error(
    "Missing function implementation: pydata");
}

void
srand(...)
{
  throw std::runtime_error(
    "Missing function implementation: srand");
}

void
zpcbdt(...)
{
  throw std::runtime_error(
    "Missing function implementation: zpcbdt");
}

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

struct common_arout
{
  int iout;

  common_arout() :
    iout(fem::int0)
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

struct common_smearz
{
  double smearp;
  double smearh;

  common_smearz() :
    smearp(fem::double0),
    smearh(fem::double0)
  {}
};

struct common_rndf77
{
  int nseed;

  common_rndf77() :
    nseed(fem::int0)
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

struct common_rndm3
{
  int iseedp;

  common_rndm3() :
    iseedp(fem::int0)
  {}
};

struct common_run
{
  int num;

  common_run() :
    num(fem::int0)
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

struct common_resdcy
{
  int nsav;
  int iksdcy;

  common_resdcy() :
    nsav(fem::int0),
    iksdcy(fem::int0)
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

struct common_cmsflag
{
  double dshadow;
  int ishadow;

  common_cmsflag() :
    dshadow(fem::double0),
    ishadow(fem::int0)
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

struct common :
  fem::common,
  common_hmain1,
  common_hparnt,
  common_ludat1,
  common_arprnt,
  common_arout,
  common_arevt,
  common_smearz,
  common_rndf77,
  common_anim,
  common_coal,
  common_snn,
  common_para2,
  common_para7,
  common_para8,
  common_rndm3,
  common_run,
  common_input1,
  common_input2,
  common_oscar1,
  common_oscar2,
  common_resdcy,
  common_phidcy,
  common_embed,
  common_cmsflag,
  common_phihj
{
  fem::cmn_sve program_ampt_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct program_ampt_save
{
  float bmax;
  float bmin;
  int i;
  int iamax;
  int ihjsed;
  int imiss;
  int ipop;
  int j;
  int k;
  int nevnt;
  int nseedr;
  fem::str<8> proj;
  fem::str<8> targ;

  program_ampt_save() :
    bmax(fem::float0),
    bmin(fem::float0),
    i(fem::int0),
    iamax(fem::int0),
    ihjsed(fem::int0),
    imiss(fem::int0),
    ipop(fem::int0),
    j(fem::int0),
    k(fem::int0),
    nevnt(fem::int0),
    nseedr(fem::int0),
    proj(fem::char0),
    targ(fem::char0)
  {}
};

//C.....driver program for A Multi-Phase Transport model
void
program_ampt(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  FEM_CMN_SVE(program_ampt);
  common_read read(cmn);
  common_write write(cmn);
  int& natt = cmn.natt;
  arr_ref<float> hipr1(cmn.hipr1, dimension(100));
  arr_ref<int> ihpr2(cmn.ihpr2, dimension(50));
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_ref<float> parj(cmn.parj, dimension(200));
  arr_ref<float> arpar1(cmn.arpar1, dimension(100));
  arr_ref<int> iaint2(cmn.iaint2, dimension(50));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  int& nseed = cmn.nseed;
  int& isoft = cmn.isoft;
  float& efrm = cmn.efrm;
  int& num = cmn.num;
  int& iap = cmn.iap;
  int& izp = cmn.izp;
  int& iat = cmn.iat;
  int& izt = cmn.izt;
  fem::str<8>& frame = cmn.frame;
  fem::str<25>& amptvn = cmn.amptvn;
  //
  float& bmax = sve.bmax;
  float& bmin = sve.bmin;
  int& i = sve.i;
  int& iamax = sve.iamax;
  int& ihjsed = sve.ihjsed;
  int& imiss = sve.imiss;
  int& ipop = sve.ipop;
  int& j = sve.j;
  int& k = sve.k;
  int& nevnt = sve.nevnt;
  int& nseedr = sve.nseedr;
  fem::str<8>& proj = sve.proj;
  fem::str<8>& targ = sve.targ;
  static const char* format_111 = "(a8)";
  static const char* format_50 =
    "(' ',/,11x,'##################################################',/,1x,10x,"
    "'#      AMPT (A Multi-Phase Transport) model      #',/,1x,10x,"
    "'#          Version ',a25,'     #',/,1x,10x,"
    "'#               10/28/2016                       #',/,1x,10x,"
    "'##################################################',/,1x,10x,' ')";
  //C
  //C     parton coalescence radii in case of string melting:
  //C     initialization value for parton cascade:
  //C     initialization value for hadron cascade:
  //Clin-4/2012-6/2009:
  //C      common/phidcy/iphidcy
  //Clin-7/2009:
  //Clin-2/2012 allow random orientation of reaction plane:
  //C
  //C****************
  cmn.io.open(24, "input.ampt")
    .status("UNKNOWN");
  cmn.io.open(12, "ana/version")
    .status("UNKNOWN");
  read(24, star), efrm;
  //C     format-read characters (for ALPHA compilers):
  read(24, format_111), frame;
  read(24, format_111), proj;
  read(24, format_111), targ;
  read(24, star), iap;
  read(24, star), izp;
  read(24, star), iat;
  read(24, star), izt;
  read(24, star), nevnt;
  read(24, star), bmin;
  read(24, star), bmax;
  //C     flag to select default AMPT or string melting:
  read(24, star), isoft;
  //C     read initialization value for hadron cascade:
  read(24, star), cmn.ntmax;
  read(24, star), cmn.dt;
  //C     parj(41) and (42) are a and b parameters in Lund string fragmentation:
  read(24, star), parj(41);
  read(24, star), parj(42);
  //C     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and
  //C     increase the net-baryon stopping in rapidity (value HIJING is 1):
  read(24, star), ipop;
  if (ipop == 1) {
    ihpr2(11) = 3;
  }
  //C     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
  read(24, star), parj(5);
  //C     shadowing flag in HIJING:
  read(24, star), ihpr2(6);
  //C     quenching flag in HIJING:
  read(24, star), ihpr2(4);
  //C     quenching rate when quenching flag is on (=1.0 GeV/fm):
  read(24, star), hipr1(14);
  //C     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV.
  read(24, star), hipr1(8);
  //C     read initialization value for parton cascade:
  read(24, star), cmn.xmu;
  read(24, star), cmn.izpc;
  read(24, star), cmn.alpha;
  //C     quark coalescence radii in momentum and space for string melting:
  read(24, star), cmn.dpcoal;
  read(24, star), cmn.drcoal;
  //C     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
  read(24, star), ihjsed;
  //C     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
  read(24, star), nseed;
  read(24, star), cmn.iseedp;
  read(24, star), cmn.iksdcy;
  read(24, star), cmn.iphidcy;
  read(24, star), cmn.ipi0dcy;
  //C     flag for OSCAR output for final partons and hadrons:
  read(24, star), cmn.ioscar;
  //Clin-5/2008     flag for perturbative treatment of deuterons:
  read(24, star), cmn.idpert;
  read(24, star), cmn.npertd;
  read(24, star), cmn.idxsec;
  //Clin-6/2009 To select events that have at least 1 high-Pt minijet parton:
  read(24, star), cmn.pttrig;
  read(24, star), cmn.maxmiss;
  read(24, star), ihpr2(2);
  read(24, star), ihpr2(5);
  //Clin-6/2009 To embed a back-to-back q/qbar pair into each event:
  read(24, star), cmn.iembed;
  read(24, star), cmn.pxqembd, cmn.pyqembd;
  read(24, star), cmn.xembd, cmn.yembd;
  read(24, star), cmn.nsembd, cmn.psembd, cmn.tmaxembd;
  //Clin-7/2009 Allow modification of nuclear shadowing:
  read(24, star), cmn.ishadow;
  read(24, star), cmn.dshadow;
  read(24, star), cmn.iphirp;
  //C
  cmn.io.close(24);
  //Clin-6/2009 ctest off turn on jet triggering:
  //C      IHPR2(3)=1
  //C     Trigger Pt of high-pt jets in HIJING:
  //C      HIPR1(10)=7.
  //C
  if (isoft == 1) {
    amptvn = "1.26t7 (Default)";
  }
  else if (isoft == 4) {
    amptvn = "2.26t7 (StringMelting)";
  }
  else {
    amptvn = "Test-Only";
  }
  write(6, format_50), amptvn;
  write(12, format_50), amptvn;
  //C     when ihjsed=11: use environment variable at run time for HIJING nseed:
  if (ihjsed == 11) {
    write(6, star), "# Read in NSEED in HIJING at run time (e.g. 20030819):";
  }
  read(6, star), nseedr;
  if (ihjsed == 11) {
    nseed = nseedr;
  }
  if (ihjsed == 11) {
    write(6, star), "#   read in: ", nseed;
    write(12, star), "# Read in NSEED in HIJING at run time:", nseed;
  }
  cmn.io.close(12);
  //Clin-5/2015 an odd number is needed for the random number generator:
  //C      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
  nseed = 2 * nseed + 1;
  //C     9/26/03 random number generator for f77 compiler:
  srand(nseed);
  //C
  //C.....turn on warning messages in nohup.out when an event is repeated:
  ihpr2(10) = 1;
  //C     string formation time:
  arpar1(1) = 0.7f;
  //C     smearp is the smearing halfwidth on parton z0,
  //C     set to 0 for now to avoid overflow in eta.
  //C     smearh is the smearing halfwidth on string production point z0.
  cmn.smearp = 0e0;
  iamax = fem::max(iap, iat);
  cmn.smearh = 1.2e0 * fem::pow(iamax, 0.3333e0) / (fem::dble(efrm) /
    2 / 0.938e0);
  cmn.nevent = nevnt;
  //C
  //C     AMPT momentum and space info at freezeout:
  cmn.io.open(16, "ana/ampt.dat")
    .status("UNKNOWN");
  cmn.io.open(14, "ana/zpc.dat")
    .status("UNKNOWN");
  //Ctest off for resonance (phi, K*) studies:
  //C      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN')
  //C      OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN')
  hijset(efrm, frame, proj, targ, iap, izp, iat, izt);
  artset();
  inizpc();
  //Clin-5/2009 ctest off:
  //C      call flowp(0)
  //C      call flowh0(NEVNT,0)
  //C      call iniflw(NEVNT,0)
  //C      call frztm(NEVNT,0)
  //C
  FEM_DO_SAFE(j, 1, nevnt) {
    iaevt = j;
    FEM_DO_SAFE(k, 1, num) {
      iarun = k;
      if (iaevt == nevnt && iarun == num) {
        cmn.iout = 1;
      }
      write(6, star), " EVENT ", j, ", RUN ", k;
      imiss = 0;
      statement_100:
      hijing(frame, bmin, bmax);
      iaint2(1) = natt;
      //C
      //Clin-6/2009 ctest off
      if (j ==  - 2) {
        write(98, star), hipr1;
        write(98, star), " ";
        write(98, star), ihpr2;
        write(98, star), " ";
        {
          write_loop wloop(cmn, 98, star);
          FEM_DO_SAFE(i, 1, 20) {
            wloop, hint1(i);
          }
        }
        write(98, star), " ";
        {
          write_loop wloop(cmn, 98, star);
          FEM_DO_SAFE(i, 21, 40) {
            wloop, hint1(i);
          }
        }
        write(98, star), " ";
        {
          write_loop wloop(cmn, 98, star);
          FEM_DO_SAFE(i, 41, 60) {
            wloop, hint1(i);
          }
        }
        write(98, star), " ";
        {
          write_loop wloop(cmn, 98, star);
          FEM_DO_SAFE(i, 61, 80) {
            wloop, hint1(i);
          }
        }
        write(98, star), " ";
        {
          write_loop wloop(cmn, 98, star);
          FEM_DO_SAFE(i, 81, 100) {
            wloop, hint1(i);
          }
        }
        write(98, star), " ";
        write(98, star), ihnt2;
      }
      //C
      //C     evaluate Npart (from primary NN collisions) for both proj and targ:
      getnp();
      //C     switch for final parton fragmentation:
      if (ihpr2(20) == 0) {
        goto statement_2000;
      }
      //C     In the unlikely case of no interaction (even after loop of 20 in HIJING),
      //C     still repeat the event to get an interaction
      //C     (this may have an additional "trigger" effect):
      if (natt == 0) {
        imiss++;
        if (imiss <= 20) {
          write(6, star), "repeated event: natt=0,j,imiss=", j, imiss;
          goto statement_100;
        }
        else {
          write(6, star), "missed event: natt=0,j=", j;
          goto statement_2000;
        }
      }
      //C.....ART initialization and run
      arini();
      arini2(k);
    }
    //C
    artan1();
    //Clin-9/2012 Analysis is not used:
    //C          CALL HJANA3
    artmn();
    //Clin-9/2012 Analysis is not used:
    //C          CALL HJANA4
    artan2();
    statement_2000:;
  }
  //C
  artout(nevnt);
  //Clin-5/2009 ctest off:
  //C       call flowh0(NEVNT,2)
  //C       call flowp(2)
  //C       call iniflw(NEVNT,2)
  //C       call frztm(NEVNT,2)
  //C
  FEM_STOP(0);
}
//C     FYI: taken file unit numbers are 10-88, 91-93;
//C     so free file unit numbers are 1-4,7-9,89,97-99.

} // namespace AMPT

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    AMPT::program_ampt);
}
