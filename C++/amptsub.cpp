#include <fem.hpp>  // Fortran EMulation library of fable module

namespace AMPT {

using namespace fem::major_types;

void distce(...) {
  throw std::runtime_error("Missing function implementation: distce");
}

void luchge(...) {
  throw std::runtime_error("Missing function implementation: luchge");
}

void rand(...) {
  throw std::runtime_error("Missing function implementation: rand");
}

void rotate(...) {
  throw std::runtime_error("Missing function implementation: rotate");
}

struct common_gg {
  float dx;
  float dy;
  float dz;
  float dpx;
  float dpy;
  float dpz;

  common_gg()
      : dx(fem::float0),
        dy(fem::float0),
        dz(fem::float0),
        dpx(fem::float0),
        dpy(fem::float0),
        dpz(fem::float0) {}
};

struct common_zz {
  int zta;
  int zpr;

  common_zz() : zta(fem::int0), zpr(fem::int0) {}
};

struct common_run {
  int num;

  common_run() : num(fem::int0) {}
};

struct common_input1 {
  int masspr;
  int massta;
  int iseed;
  int iavoid;
  float dt;

  common_input1()
      : masspr(fem::int0),
        massta(fem::int0),
        iseed(fem::int0),
        iavoid(fem::int0),
        dt(fem::float0) {}
};

struct common_input2 {
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

  common_input2()
      : ilab(fem::int0),
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
        kmul(fem::int0) {}
};

struct common_input3 {
  float plab;
  float elab;
  float zeropt;
  float b0;
  float bi;
  float bm;
  float dencut;
  float cycbox;

  common_input3()
      : plab(fem::float0),
        elab(fem::float0),
        zeropt(fem::float0),
        b0(fem::float0),
        bi(fem::float0),
        bm(fem::float0),
        dencut(fem::float0),
        cycbox(fem::float0) {}
};

struct common_imulst {
  int iperts;

  common_imulst() : iperts(fem::int0) {}
};

struct common_coal {
  double dpcoal;
  double drcoal;
  double ecritl;

  common_coal()
      : dpcoal(fem::double0), drcoal(fem::double0), ecritl(fem::double0) {}
};

struct common_anim {
  int nevent;
  int isoft;
  int isflag;
  int izpc;

  common_anim()
      : nevent(fem::int0),
        isoft(fem::int0),
        isflag(fem::int0),
        izpc(fem::int0) {}
};

struct common_para7 {
  int ioscar;
  int nsmbbbar;
  int nsmmeson;

  common_para7()
      : ioscar(fem::int0), nsmbbbar(fem::int0), nsmmeson(fem::int0) {}
};

struct common_embed {
  int iembed;
  int nsembd;
  float pxqembd;
  float pyqembd;
  float xembd;
  float yembd;
  float psembd;
  float tmaxembd;
  float phidecomp;

  common_embed()
      : iembed(fem::int0),
        nsembd(fem::int0),
        pxqembd(fem::float0),
        pyqembd(fem::float0),
        xembd(fem::float0),
        yembd(fem::float0),
        psembd(fem::float0),
        tmaxembd(fem::float0),
        phidecomp(fem::float0) {}
};

struct common_xyembed {
  static const int nxymax = 10001;

  int nxyjet;
  arr<float, 2> xyjet;

  common_xyembed()
      : nxyjet(fem::int0), xyjet(dimension(nxymax, 2), fem::fill0) {}
};

const int common_xyembed::nxymax;

struct common_arprnt {
  arr<float> arpar1;
  arr<int> iapar2;
  arr<float> arint1;
  arr<int> iaint2;

  common_arprnt()
      : arpar1(dimension(100), fem::fill0),
        iapar2(dimension(50), fem::fill0),
        arint1(dimension(100), fem::fill0),
        iaint2(dimension(50), fem::fill0) {}
};

struct common_arprc {
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

  common_arprc()
      : itypar(dimension(maxstr), fem::fill0),
        gxar(dimension(maxstr), fem::fill0),
        gyar(dimension(maxstr), fem::fill0),
        gzar(dimension(maxstr), fem::fill0),
        ftar(dimension(maxstr), fem::fill0),
        pxar(dimension(maxstr), fem::fill0),
        pyar(dimension(maxstr), fem::fill0),
        pzar(dimension(maxstr), fem::fill0),
        pear(dimension(maxstr), fem::fill0),
        xmar(dimension(maxstr), fem::fill0) {}
};

const int common_arprc::maxstr;

struct common_dpert {
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

  common_dpert()
      : dpertt(dimension(maxstr, maxr), fem::fill0),
        dpertp(dimension(maxstr), fem::fill0),
        dplast(dimension(maxstr), fem::fill0),
        dpdcy(dimension(maxstr), fem::fill0),
        dpdpi(dimension(maxstr, maxr), fem::fill0),
        dpt(dimension(maxstr, maxr), fem::fill0),
        dpp1(dimension(maxstr, maxr), fem::fill0),
        dppion(dimension(maxstr, maxr), fem::fill0) {}
};

const int common_dpert::maxstr;
const int common_dpert::maxr;

struct common_smearz {
  double smearp;
  double smearh;

  common_smearz() : smearp(fem::double0), smearh(fem::double0) {}
};

struct common_rndf77 {
  int nseed;

  common_rndf77() : nseed(fem::int0) {}
};

struct common_para8 {
  int idpert;
  int npertd;
  int idxsec;

  common_para8() : idpert(fem::int0), npertd(fem::int0), idxsec(fem::int0) {}
};

struct common_nzpc {
  int nattzp;

  common_nzpc() : nattzp(fem::int0) {}
};

struct common_hparnt {
  arr<float> hipr1;
  arr<int> ihpr2;
  arr<float> hint1;
  arr<int> ihnt2;

  common_hparnt()
      : hipr1(dimension(100), fem::fill0),
        ihpr2(dimension(50), fem::fill0),
        hint1(dimension(100), fem::fill0),
        ihnt2(dimension(50), fem::fill0) {}
};

struct common_arerc1 {
  static const int maxr = 1;

  arr<int> multi1;

  common_arerc1() : multi1(dimension(maxr), fem::fill0) {}
};

const int common_arerc1::maxr;

struct common_arprc1 {
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

  common_arprc1()
      : ityp1(dimension(maxstr, maxr), fem::fill0),
        gx1(dimension(maxstr, maxr), fem::fill0),
        gy1(dimension(maxstr, maxr), fem::fill0),
        gz1(dimension(maxstr, maxr), fem::fill0),
        ft1(dimension(maxstr, maxr), fem::fill0),
        px1(dimension(maxstr, maxr), fem::fill0),
        py1(dimension(maxstr, maxr), fem::fill0),
        pz1(dimension(maxstr, maxr), fem::fill0),
        ee1(dimension(maxstr, maxr), fem::fill0),
        xm1(dimension(maxstr, maxr), fem::fill0) {}
};

const int common_arprc1::maxstr;
const int common_arprc1::maxr;

struct common_tdecay {
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float> tfdcy;
  arr<float, 2> tfdpi;
  arr<float> tft;

  common_tdecay()
      : tfdcy(dimension(maxstr), fem::fill0),
        tfdpi(dimension(maxstr, maxr), fem::fill0),
        tft(dimension(maxstr), fem::fill0) {}
};

const int common_tdecay::maxstr;
const int common_tdecay::maxr;

struct common_aa {
  static const int maxstr = 150001;

  arr<float, 2> r;

  common_aa() : r(dimension(3, maxstr), fem::fill0) {}
};

const int common_aa::maxstr;

struct common_bb {
  static const int maxstr = 150001;

  arr<float, 2> p;

  common_bb() : p(dimension(3, maxstr), fem::fill0) {}
};

const int common_bb::maxstr;

struct common_cc {
  static const int maxstr = 150001;

  arr<float> e;

  common_cc() : e(dimension(maxstr), fem::fill0) {}
};

const int common_cc::maxstr;

struct common_ee {
  static const int maxstr = 150001;

  arr<int> id;
  arr<int> lb;

  common_ee()
      : id(dimension(maxstr), fem::fill0), lb(dimension(maxstr), fem::fill0) {}
};

const int common_ee::maxstr;

struct common_bg {
  float betax;
  float betay;
  float betaz;
  float gamma;

  common_bg()
      : betax(fem::float0),
        betay(fem::float0),
        betaz(fem::float0),
        gamma(fem::float0) {}
};

struct common_nn {
  int nnn;

  common_nn() : nnn(fem::int0) {}
};

struct common_pa {
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 3> rpion;

  common_pa() : rpion(dimension(3, maxstr, maxr), fem::fill0) {}
};

const int common_pa::maxstr;
const int common_pa::maxr;

struct common_pb {
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 3> ppion;

  common_pb() : ppion(dimension(3, maxstr, maxr), fem::fill0) {}
};

const int common_pb::maxstr;
const int common_pb::maxr;

struct common_pc {
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<float, 2> epion;

  common_pc() : epion(dimension(maxstr, maxr), fem::fill0) {}
};

const int common_pc::maxstr;
const int common_pc::maxr;

struct common_pd {
  static const int maxstr = 150001;
  static const int maxr = 1;

  arr<int, 2> lpion;

  common_pd() : lpion(dimension(maxstr, maxr), fem::fill0) {}
};

const int common_pd::maxstr;
const int common_pd::maxr;

struct common_arana1 {
  arr<float> dy1ntb;
  arr<float> dy1ntp;
  arr<float> dy1hm;
  arr<float> dy1kp;
  arr<float> dy1km;
  arr<float> dy1k0s;
  arr<float> dy1la;
  arr<float> dy1lb;
  arr<float> dy1phi;
  arr<float> dm1pip;
  arr<float> dm1pim;
  arr<float> dmt1pr;
  arr<float> dmt1pb;
  arr<float> dmt1kp;
  arr<float> dm1km;
  arr<float> dm1k0s;
  arr<float> dmt1la;
  arr<float> dmt1lb;
  arr<float> dy1msn;
  arr<float> dy1pip;
  arr<float> dy1pim;
  arr<float> dy1pi0;
  arr<float> dy1pr;
  arr<float> dy1pb;
  arr<float> dy1neg;
  arr<float> dy1ch;
  arr<float> de1neg;
  arr<float> de1ch;

  common_arana1()
      : dy1ntb(dimension(50), fem::fill0),
        dy1ntp(dimension(50), fem::fill0),
        dy1hm(dimension(50), fem::fill0),
        dy1kp(dimension(50), fem::fill0),
        dy1km(dimension(50), fem::fill0),
        dy1k0s(dimension(50), fem::fill0),
        dy1la(dimension(50), fem::fill0),
        dy1lb(dimension(50), fem::fill0),
        dy1phi(dimension(50), fem::fill0),
        dm1pip(dimension(50), fem::fill0),
        dm1pim(dimension(50), fem::fill0),
        dmt1pr(dimension(50), fem::fill0),
        dmt1pb(dimension(50), fem::fill0),
        dmt1kp(dimension(50), fem::fill0),
        dm1km(dimension(50), fem::fill0),
        dm1k0s(dimension(50), fem::fill0),
        dmt1la(dimension(50), fem::fill0),
        dmt1lb(dimension(50), fem::fill0),
        dy1msn(dimension(50), fem::fill0),
        dy1pip(dimension(50), fem::fill0),
        dy1pim(dimension(50), fem::fill0),
        dy1pi0(dimension(50), fem::fill0),
        dy1pr(dimension(50), fem::fill0),
        dy1pb(dimension(50), fem::fill0),
        dy1neg(dimension(50), fem::fill0),
        dy1ch(dimension(50), fem::fill0),
        de1neg(dimension(50), fem::fill0),
        de1ch(dimension(50), fem::fill0) {}
};

struct common_arana2 {
  arr<float> dy2ntb;
  arr<float> dy2ntp;
  arr<float> dy2hm;
  arr<float> dy2kp;
  arr<float> dy2km;
  arr<float> dy2k0s;
  arr<float> dy2la;
  arr<float> dy2lb;
  arr<float> dy2phi;
  arr<float> dm2pip;
  arr<float> dm2pim;
  arr<float> dmt2pr;
  arr<float> dmt2pb;
  arr<float> dmt2kp;
  arr<float> dm2km;
  arr<float> dm2k0s;
  arr<float> dmt2la;
  arr<float> dmt2lb;
  arr<float> dy2msn;
  arr<float> dy2pip;
  arr<float> dy2pim;
  arr<float> dy2pi0;
  arr<float> dy2pr;
  arr<float> dy2pb;
  arr<float> dy2neg;
  arr<float> dy2ch;
  arr<float> de2neg;
  arr<float> de2ch;

  common_arana2()
      : dy2ntb(dimension(50), fem::fill0),
        dy2ntp(dimension(50), fem::fill0),
        dy2hm(dimension(50), fem::fill0),
        dy2kp(dimension(50), fem::fill0),
        dy2km(dimension(50), fem::fill0),
        dy2k0s(dimension(50), fem::fill0),
        dy2la(dimension(50), fem::fill0),
        dy2lb(dimension(50), fem::fill0),
        dy2phi(dimension(50), fem::fill0),
        dm2pip(dimension(50), fem::fill0),
        dm2pim(dimension(50), fem::fill0),
        dmt2pr(dimension(50), fem::fill0),
        dmt2pb(dimension(50), fem::fill0),
        dmt2kp(dimension(50), fem::fill0),
        dm2km(dimension(50), fem::fill0),
        dm2k0s(dimension(50), fem::fill0),
        dmt2la(dimension(50), fem::fill0),
        dmt2lb(dimension(50), fem::fill0),
        dy2msn(dimension(50), fem::fill0),
        dy2pip(dimension(50), fem::fill0),
        dy2pim(dimension(50), fem::fill0),
        dy2pi0(dimension(50), fem::fill0),
        dy2pr(dimension(50), fem::fill0),
        dy2pb(dimension(50), fem::fill0),
        dy2neg(dimension(50), fem::fill0),
        dy2ch(dimension(50), fem::fill0),
        de2neg(dimension(50), fem::fill0),
        de2ch(dimension(50), fem::fill0) {}
};

struct common_para1 {
  int mul;

  common_para1() : mul(fem::int0) {}
};

struct common_hjcrdn {
  arr<float, 2> yp;
  arr<float, 2> yt;

  common_hjcrdn()
      : yp(dimension(3, 300), fem::fill0), yt(dimension(3, 300), fem::fill0) {}
};

struct common_hjjet1 {
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

  common_hjjet1()
      : npj(dimension(300), fem::fill0),
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
        pjtm(dimension(300, 500), fem::fill0) {}
};

struct common_hjjet2 {
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

  common_hjjet2()
      : nsg(fem::int0),
        njsg(dimension(maxstr), fem::fill0),
        iasg(dimension(maxstr, 3), fem::fill0),
        k1sg(dimension(maxstr, 100), fem::fill0),
        k2sg(dimension(maxstr, 100), fem::fill0),
        pxsg(dimension(maxstr, 100), fem::fill0),
        pysg(dimension(maxstr, 100), fem::fill0),
        pzsg(dimension(maxstr, 100), fem::fill0),
        pesg(dimension(maxstr, 100), fem::fill0),
        pmsg(dimension(maxstr, 100), fem::fill0) {}
};

const int common_hjjet2::maxstr;

struct common_prec1 {
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

  common_prec1()
      : gx0(dimension(maxptn), fem::fill0),
        gy0(dimension(maxptn), fem::fill0),
        gz0(dimension(maxptn), fem::fill0),
        ft0(dimension(maxptn), fem::fill0),
        px0(dimension(maxptn), fem::fill0),
        py0(dimension(maxptn), fem::fill0),
        pz0(dimension(maxptn), fem::fill0),
        e0(dimension(maxptn), fem::fill0),
        xmass0(dimension(maxptn), fem::fill0),
        ityp0(dimension(maxptn), fem::fill0) {}
};

const int common_prec1::maxptn;

struct common_arevt {
  int iaevt;
  int iarun;
  int miss;

  common_arevt() : iaevt(fem::int0), iarun(fem::int0), miss(fem::int0) {}
};

struct common_arout {
  int iout;

  common_arout() : iout(fem::int0) {}
};

struct common_prec2 {
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

  common_prec2()
      : gx5(dimension(maxptn), fem::fill0),
        gy5(dimension(maxptn), fem::fill0),
        gz5(dimension(maxptn), fem::fill0),
        ft5(dimension(maxptn), fem::fill0),
        px5(dimension(maxptn), fem::fill0),
        py5(dimension(maxptn), fem::fill0),
        pz5(dimension(maxptn), fem::fill0),
        e5(dimension(maxptn), fem::fill0),
        xmass5(dimension(maxptn), fem::fill0),
        ityp5(dimension(maxptn), fem::fill0) {}
};

const int common_prec2::maxptn;

struct common_ilist8 {
  static const int maxptn = 400001;

  arr<int> lstrg1;
  arr<int> lpart1;

  common_ilist8()
      : lstrg1(dimension(maxptn), fem::fill0),
        lpart1(dimension(maxptn), fem::fill0) {}
};

const int common_ilist8::maxptn;

struct common_srec1 {
  int nsp;
  int nst;
  int nsi;

  common_srec1() : nsp(fem::int0), nst(fem::int0), nsi(fem::int0) {}
};

struct common_srec2 {
  static const int maxstr = 150001;

  arr<double> ataui;
  arr<double> zt1;
  arr<double> zt2;
  arr<double> zt3;

  common_srec2()
      : ataui(dimension(maxstr), fem::fill0),
        zt1(dimension(maxstr), fem::fill0),
        zt2(dimension(maxstr), fem::fill0),
        zt3(dimension(maxstr), fem::fill0) {}
};

const int common_srec2::maxstr;

struct common_soft {
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

  common_soft()
      : pxsgs(dimension(maxstr, 3), fem::fill0),
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
        njsgs(dimension(maxstr), fem::fill0) {}
};

const int common_soft::maxstr;

struct common_iflow {
  double v2i;
  double eti;
  double xmulti;
  double v2mi;
  double s2mi;
  double xmmult;
  double v2bi;
  double s2bi;
  double xbmult;

  common_iflow()
      : v2i(fem::double0),
        eti(fem::double0),
        xmulti(fem::double0),
        v2mi(fem::double0),
        s2mi(fem::double0),
        xmmult(fem::double0),
        v2bi(fem::double0),
        s2bi(fem::double0),
        xbmult(fem::double0) {}
};

struct common_fflow {
  float v2f;
  float etf;
  float xmultf;
  float v2fpi;
  float xmulpi;

  common_fflow()
      : v2f(fem::float0),
        etf(fem::float0),
        xmultf(fem::float0),
        v2fpi(fem::float0),
        xmulpi(fem::float0) {}
};

struct common_strg {
  static const int maxstr = 150001;

  arr<int> np;

  common_strg() : np(dimension(maxstr), fem::fill0) {}
};

const int common_strg::maxstr;

struct common_frzprc {
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

  common_frzprc()
      : gxfrz(dimension(maxptn), fem::fill0),
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
        itlast(fem::int0) {}
};

const int common_frzprc::maxptn;

struct common : fem::common,
                common_gg,
                common_zz,
                common_run,
                common_input1,
                common_input2,
                common_input3,
                common_imulst,
                common_coal,
                common_anim,
                common_para7,
                common_embed,
                common_xyembed,
                common_arprnt,
                common_arprc,
                common_dpert,
                common_smearz,
                common_rndf77,
                common_para8,
                common_nzpc,
                common_hparnt,
                common_arerc1,
                common_arprc1,
                common_tdecay,
                common_aa,
                common_bb,
                common_cc,
                common_ee,
                common_bg,
                common_nn,
                common_pa,
                common_pb,
                common_pc,
                common_pd,
                common_arana1,
                common_arana2,
                common_para1,
                common_hjcrdn,
                common_hjjet1,
                common_hjjet2,
                common_prec1,
                common_arevt,
                common_arout,
                common_prec2,
                common_ilist8,
                common_srec1,
                common_srec2,
                common_soft,
                common_iflow,
                common_fflow,
                common_strg,
                common_frzprc {
  fem::cmn_sve artset_sve;
  fem::cmn_sve addhad_sve;
  fem::cmn_sve arini1_sve;
  fem::cmn_sve arindx_sve;
  fem::cmn_sve artord_sve;
  fem::cmn_sve arini_sve;
  fem::cmn_sve arini2_sve;
  fem::cmn_sve iarflv_sve;
  fem::cmn_sve blockdata_ardata_sve;
  fem::cmn_sve x2kaon_sve;
  fem::cmn_sve pinsg0_sve;
  fem::cmn_sve aknel_sve;
  fem::cmn_sve akpel_sve;
  fem::cmn_sve aknsgm_sve;
  fem::cmn_sve akpsgm_sve;
  fem::cmn_sve akplam_sve;
  fem::cmn_sve fstate_sve;
  fem::cmn_sve nnkaon_sve;
  fem::cmn_sve lorntz_sve;
  fem::cmn_sve npik_sve;
  fem::cmn_sve pihypn_sve;
  fem::cmn_sve kaonn_sve;
  fem::cmn_sve newka_sve;
  fem::cmn_sve aknpsg_sve;
  fem::cmn_sve artan1_sve;
  fem::cmn_sve artan2_sve;
  fem::cmn_sve artout_sve;
  fem::cmn_sve hjana1_sve;
  fem::cmn_sve hjan1b_sve;
  fem::cmn_sve hjan1a_sve;
  fem::cmn_sve hjan2a_sve;
  fem::cmn_sve hjan2b_sve;
  fem::cmn_sve hjana2_sve;
  fem::cmn_sve hjana3_sve;
  fem::cmn_sve hjana4_sve;
  fem::cmn_sve zpstrg_sve;

  common(int argc, char const* argv[]) : fem::common(argc, argv) {}
};

struct artset_save {
  int iplab;
  int ixy;

  artset_save() : iplab(fem::int0), ixy(fem::int0) {}
};

// C....................amptsub.f
// C.....this file contains 4 sections:
// C.....1. ART subroutines;
// C.....2. ART functions;
// C.....3. ART block data;
// C.....4. subprocesses borrowed from other codes.
// C.....5. the previous artana.f
// C.....6. the previous zpcsub.f
// C.....7. subroutine getnp
// C.....Note that Parts1-4 are the previous artsub.f
// C
// C=======================================================================
// C.....subroutine to set up ART parameters and analysis files
// C.....before looping different events
void artset(common& cmn) {
  FEM_CMN_SVE(artset);
  common_read read(cmn);
  common_write write(cmn);
  // COMMON input2
  int& icoll = cmn.icoll;
  int& ipot = cmn.ipot;
  int& imomen = cmn.imomen;
  // COMMON input3
  float& plab = cmn.plab;
  float& elab = cmn.elab;
  // COMMON anim
  int& isoft = cmn.isoft;
  // COMMON para7
  int& ioscar = cmn.ioscar;
  // COMMON embed
  int& iembed = cmn.iembed;
  // COMMON xyembed
  int& nxyjet = cmn.nxyjet;
  const int nxymax = 10001;
  arr_ref<float, 2> xyjet(cmn.xyjet, dimension(nxymax, 2));
  //
  // SAVE
  int& iplab = sve.iplab;
  int& ixy = sve.ixy;
  //
  // C
  // Clin-10/03/03
  // C     "SAVE   " (without argument) is used for most subroutines and
  // functions, C     this is important for the success when using "f77" to
  // compile: Cc      SAVE /gg/ Cc      SAVE /zz/ Cc      SAVE /RUN/ Cc      SAVE
  // /input1/ Cc      SAVE /INPUT2/ Cc      SAVE /INPUT3/ Cc      SAVE /imulst/
  // Clin-10/03/03  ecritl: local energy density below which a parton
  // C     will freeze out (in GeV/fm^3), for improvements on string melting,
  // C     not used in this version of AMPT:
  // Clin-4/2008
  // C      data ecritl/1.d0/
  cmn.ecritl = 1.e0;
  // C
  // C     combine ART initialization into ampt.ini:
  // C     (Note that the following values are relics from the old ART
  // structure) C.....input parameter file C      OPEN(13, FILE = 'art1.ini',
  // STATUS = 'UNKNOWN') C      READ (13, *) MASSTA, ZTA
  cmn.massta = 1;
  cmn.zta = 1;
  // C      write(12,*) massta, zta, ' massta, zta'
  // C      READ (13, *) MASSPR, ZPR
  cmn.masspr = 1;
  cmn.zpr = 1;
  // C      write(12,*) masspr, zpr, ' masspr, zpr'
  // C      READ (13, *) PLAB, IPLAB
  plab = 14.6f;
  iplab = 2;
  // C      write(12,*) plab, iplab, ' plab, iplab'
  const float amu = 0.9383f;
  if (iplab == 2) {
    elab = fem::sqrt(fem::pow2(plab) + fem::pow2(amu)) - amu;
  } else {
    elab = plab;
  }
  elab = elab * 1000.f;
  // C      READ (13, *) ZEROPT
  cmn.zeropt = 0.f;
  // C      write(12,*) zeropt, ' zeropt'
  // Clin-10/03/03 ISEED was used as a seed for random number inside ART,
  // C     not used in AMPT:
  cmn.iseed = 700721;
  // C     0/1: (Normal or Perturbative) multistrange partice production.
  // C     Perturbative option is disabled for now:
  cmn.iperts = 0;
  // C      READ (13, *) MANYB, B0, BI, BM
  // C     2/04/00 MANYB MUST BE SET TO 1 !
  // C     in order to skip impact parameter setting by ART, then B0 has no
  // effect.
  cmn.manyb = 1;
  cmn.b0 = 1;
  cmn.bi = 0;
  cmn.bm = 0;
  // C      write(12,*) manyb, b0, bi, bm, ' manyb, b0, bi, bm'
  // C      READ (13, *) ISEED
  // C      write(12,*) iseed, ' iseed'
  // C      READ (13, *) DT
  // C      write(12,*) dt, ' dt'
  // C      READ (13, *) NTMAX
  // C      write(12,*) ntmax, ' ntmax'
  // C      READ (13, *) ICOLL
  icoll = -1;
  // C      write(12,*) icoll, ' icoll'
  // C      READ (13, *) NUM
  // C     2/11/03 run events without test particles for now:
  cmn.num = 1;
  // C      write(12,*) num, ' num'
  // C      READ (13, *) INSYS
  cmn.insys = 1;
  // C      write(12,*) insys, ' insys'
  // C      READ (13, *) IPOT
  ipot = 3;
  // C      write(12,*) ipot, ' ipot'
  // C      READ (13, *) MODE
  cmn.mode = 0;
  if (icoll == -1) {
    ipot = 0;
  }
  // C      write(12,*) mode, ' mode'
  // C      READ (13, *) DX, DY, DZ
  cmn.dx = 2.73f;
  cmn.dy = 2.73f;
  cmn.dz = 2.73f;
  // C      write(12,*) dx,dy,dz,' dx,dy,dz'
  // C      READ (13, *) DPX, DPY, DPZ
  cmn.dpx = 0.6f;
  cmn.dpy = 0.6f;
  cmn.dpz = 0.6f;
  // C      write(12,*) dpx,dpy,dpz,' dpx,dpy,dpz'
  // C      READ (13, *) IAVOID
  cmn.iavoid = 1;
  // C      write(12,*) iavoid, ' iavoid'
  // C      READ (13, *) IMOMEN
  imomen = 1;
  // C      write(12,*) imomen, ' imomen'
  if (icoll == -1) {
    imomen = 3;
  }
  // C      READ (13, *) NFREQ
  cmn.nfreq = 10;
  // C      write(12,*) nfreq, ' nfreq'
  // C      READ (13, *) ICFLOW
  cmn.icflow = 0;
  // C      write(12,*) ICFLOW, ' ICFLOW'
  // C      READ (13, *) ICRHO
  cmn.icrho = 0;
  // C      write(12,*) ICRHO, ' ICRHO'
  // C      READ (13, *) ICOU
  cmn.icou = 0;
  // C      write(12,*)icou, ' icou'
  // C kaon potential control parameter
  // C KMUL IS A MULTIPLIER TO THE STANDARD K-N SCATTERING LENGTH
  // C      READ (13, *) KPOTEN, KMUL
  cmn.kpoten = 0;
  cmn.kmul = 1;
  // C      write(12,*)kpoten,kmul, ' kpoten, kmul'
  // C mean field control parameter FOR BARYONS
  // C no mean filed is used for baryons if their
  // C local density is higher than dencut.
  // C      READ (13, *) DENCUT
  cmn.dencut = 15;
  // C      write(12,*)dencut, ' dencut'
  // C test reactions in a box of side-length cycbox
  // C input cycbox
  // C      READ (13, *) CYCBOX
  cmn.cycbox = 0;
  // C      write(12,*) cycbox, ' cycbox'
  // C
  // Clin-5b/2008
  // C      if(ioscar.eq.2) then
  if (ioscar == 2 || ioscar == 3) {
    cmn.io.open(92, "ana/parton-initial-afterPropagation.dat")
        .status("UNKNOWN");
  }

  if (ioscar == 3) {
    // Clin-6/2009 write out full parton collision history:
    cmn.io.open(95, "ana/parton-collisionsHistory.dat").status("UNKNOWN");
    // Clin-6/2009 write out initial minijet information:
    cmn.io.open(96, "ana/minijet-initial-beforePropagation.dat")
        .status("UNKNOWN");
    // Clin-6/2009 write out parton info after coalescence:
    if (isoft == 4 || isoft == 5) {
      cmn.io.open(85, "ana/parton-after-coalescence.dat").status("UNKNOWN");
    }
  }
  // Clin-6/2009 write out initial transverse positions of initial nucleons:
  cmn.io.open(94, "ana/npart-xy.dat").status("UNKNOWN");
  // C
  // Clin-8/2009 In case that random positions are used to embed high-Pt jets:
  if (iembed == 3 || iembed == 4) {
    cmn.io.open(97, "embed-jet-xy.txt").status("UNKNOWN");
    read(97, star), nxyjet;
    // C     Save positions in array to reuse when embedding more jet pairs
    // C     than the number of entries in the position file:
    if (cmn.nevent > nxyjet) {
      if (nxyjet > nxymax) {
        write(6, star),
            "Too many lines in embed-jet-xy.txt:  increase value of the "
            "paramete"
            "r nxymax";
        FEM_STOP(0);
      } else if (nxyjet <= 0) {
        write(6, star), "Check number of entries in embed-jet-xy.txt";
        FEM_STOP(0);
      }
      FEM_DO_SAFE(ixy, 1, nxyjet) {
        read(97, star), xyjet(ixy, 1), xyjet(ixy, 2);
      }
    }
  }
  // C
}

// C
// Clin-10/01/03 random number generator for f77:
float ranart(int const& /* nseed */) {
  float return_value = fem::float0;
  // Clin-4/2008 ran(nseed) is renamed to avoid conflict with system functions:
  // C      ran=rand()
  return_value = rand(0);
  // C     one may also use the following random number generator in
  // PYTHIA/JETSET: C      ranart=rlu(0)
  return return_value;
}

// C
// Clin-8/2014 define function asinh():
float asinh(float const& x) {
  float return_value = fem::float0;
  if (x > 0) {
    return_value = fem::alog(x + fem::sqrt(fem::pow2(x) + 1.f));
  } else {
    // C     a la suggestion de YP Liu:
    return_value = -fem::alog(-x + fem::sqrt(fem::pow2(x) + 1.f));
  }
  return return_value;
}

struct addhad_save {
  int i;
  int nadd;
  int np0;
  float rap;
  float tau0;
  float taui;
  float vx;
  float vy;
  float zsmear;

  addhad_save()
      : i(fem::int0),
        nadd(fem::int0),
        np0(fem::int0),
        rap(fem::float0),
        tau0(fem::float0),
        taui(fem::float0),
        vx(fem::float0),
        vy(fem::float0),
        zsmear(fem::float0) {}
};

// C
// Clin-3/2009
// C     Initialize hadron weights;
// C     Can add initial hadrons before the hadron cascade starts (but after
// ZPC).
void addhad(common& cmn) {
  FEM_CMN_SVE(addhad);
  common_write write(cmn);
  // COMMON arprnt
  arr_cref<float> arpar1(cmn.arpar1, dimension(100));
  arr_ref<int> iaint2(cmn.iaint2, dimension(50));
  // COMMON arprc
  const int maxstr = 150001;
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
  // COMMON dpert
  arr_ref<float> dpertp(cmn.dpertp, dimension(maxstr));
  // COMMON rndf77
  int& nseed = cmn.nseed;
  // COMMON para8
  int& idpert = cmn.idpert;
  //
  // SAVE
  int& i = sve.i;
  int& nadd = sve.nadd;
  int& np0 = sve.np0;
  float& rap = sve.rap;
  float& tau0 = sve.tau0;
  float& taui = sve.taui;
  float& vx = sve.vx;
  float& vy = sve.vy;
  float& zsmear = sve.zsmear;
  //
  // C
  // C     All hadrons at the start of hadron cascade have the weight of 1
  // C     except those inserted by the user in this subroutine:
  np0 = iaint2(1);
  FEM_DO_SAFE(i, 1, np0) { dpertp(i) = 1.f; }
  // C     Specify number, species, weight, initial x,p,m for inserted hadrons
  // here:
  nadd = 0;
  tau0 = arpar1(1);
  const float xmd = 1.8756f;
  FEM_DO_SAFE(i, np0 + 1, np0 + nadd) {
    itypar(i) = 42;
    // Clin-5/2012 fix type mismatch:
    // C         dpertp(I)=1d0/dble(nadd)
    dpertp(i) = 1.f / fem::ffloat(nadd);
    gxar(i) = 5.f * (1.f - 2.f * ranart(nseed));
    gyar(i) = 5.f * (1.f - 2.f * ranart(nseed));
    gzar(i) = 2.f * (1.f - 2.f * ranart(nseed));
    ftar(i) = 0.f;
    pxar(i) = 1.f;
    pyar(i) = 0.f;
    pzar(i) = 1.f;
    xmar(i) = xmd;
    // C
    pear(i) = fem::sqrt(fem::pow2(pxar(i)) + fem::pow2(pyar(i)) +
                        fem::pow2(pzar(i)) + fem::pow2(xmar(i)));
    // Clin-9/2012 determine rapidity more generally:
    // C         RAP=0.5*alog((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
    rap = asinh(pzar(i) / fem::sqrt(fem::pow2(xmar(i)) + fem::pow2(pxar(i)) +
                                    fem::pow2(pyar(i))));
    // C
    vx = pxar(i) / pear(i);
    vy = pyar(i) / pear(i);
    // C.....give initial formation time shift and boost according to rapidity:
    taui = ftar(i) + tau0;
    ftar(i) = taui * fem::cosh(rap);
    gxar(i) += vx * tau0 * fem::cosh(rap);
    gyar(i) += vy * tau0 * fem::cosh(rap);
    // C     Allow the intial z-position to be different from the Bjorken
    // picture:
    gzar(i) += taui * fem::sinh(rap);
    // C         GZAR(I)=TAUI*SINH(RAP)
    zsmear = fem::sngl(cmn.smearh) * (2.f * ranart(nseed) - 1.f);
    gzar(i) += zsmear;
  }
  iaint2(1) += nadd;
  // C
  if (nadd >= 1 && idpert != 1 && idpert != 2) {
    write(16, star),
        "IDPERT must be 1 or 2 to add initial hadrons, set NPERTD to 0 if you "
        "do"
        " not need perturbative deuterons";
    FEM_STOP(0);
  }
  if (iaint2(1) > maxstr) {
    write(16, star), "Too many initial hadrons, array size is exceeded!";
    FEM_STOP(0);
  }
  // C
}

struct arini1_save {
  int i;
  int np;
  float rap;
  float tau0;
  float taui;
  float vx;
  float vy;
  float zsmear;

  arini1_save()
      : i(fem::int0),
        np(fem::int0),
        rap(fem::float0),
        tau0(fem::float0),
        taui(fem::float0),
        vx(fem::float0),
        vy(fem::float0),
        zsmear(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....subroutine to generate formation time and position at formation time
// C.....from read-in initial conditions with an averaged formation proper
// C.....time.
// C
void arini1(common& cmn) {
  FEM_CMN_SVE(arini1);
  common_write write(cmn);
  // COMMON arprnt
  arr_cref<float> arpar1(cmn.arpar1, dimension(100));
  arr_cref<int> iaint2(cmn.iaint2, dimension(50));
  // COMMON arprc
  const int maxstr = 150001;
  arr_cref<int> itypar(cmn.itypar, dimension(maxstr));
  arr_ref<float> gxar(cmn.gxar, dimension(maxstr));
  arr_ref<float> gyar(cmn.gyar, dimension(maxstr));
  arr_ref<float> gzar(cmn.gzar, dimension(maxstr));
  arr_ref<float> ftar(cmn.ftar, dimension(maxstr));
  arr_cref<float> pxar(cmn.pxar, dimension(maxstr));
  arr_cref<float> pyar(cmn.pyar, dimension(maxstr));
  arr_cref<float> pzar(cmn.pzar, dimension(maxstr));
  arr_cref<float> pear(cmn.pear, dimension(maxstr));
  arr_cref<float> xmar(cmn.xmar, dimension(maxstr));
  // COMMON anim
  int& isoft = cmn.isoft;
  // COMMON nzpc
  int& nattzp = cmn.nattzp;
  // COMMON hparnt
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  // COMMON para8
  int& idpert = cmn.idpert;
  //
  // SAVE
  int& i = sve.i;
  int& np = sve.np;
  float& rap = sve.rap;
  float& tau0 = sve.tau0;
  float& taui = sve.taui;
  float& vx = sve.vx;
  float& vy = sve.vy;
  float& zsmear = sve.zsmear;
  //
  // C
  // C.....before invoking ARINI1:
  // C.....ARPAR1(1), IAINT2(1) must be set:
  // C
  // Cc      SAVE /ARPRNT/
  // Cc      SAVE /ARPRC/
  // Cc      SAVE /smearz/
  // Cc      SAVE /input1/
  // Cc      SAVE /anim/
  // Cc      SAVE /nzpc/
  // Cc      SAVE /HPARNT/
  // Cc      SAVE /RNDF77/
  // C
  // Clin-5/2008 for perturbatively-produced hadrons (currently only deuterons):
  cmn.io.open(91, "ana/deuteron_processes.dat").status("UNKNOWN");
  if (idpert == 1 || idpert == 2) {
    cmn.io.open(90, "ana/ampt_pert.dat").status("UNKNOWN");
  }
  // C.....generate formation time and position at formation time.
  tau0 = arpar1(1);
  np = iaint2(1);
  // Clin-7/10/01     initial positions already given for hadrons
  // C     formed from partons inside ZPC (from string melting):
  if (isoft == 3 || isoft == 4 || isoft == 5) {
    // Clin-8/2015 fixed a bug that may skip "dpertp(I)=1." in addhad and
    // C     cause the first few events to be missing in ampt.dat
    // C     (mostly for low-multiplicity events such as PP collisions):
    // C         if(NP.le.nattzp) return
    if (np > nattzp) {
      // C
      FEM_DO_SAFE(i, nattzp + 1, np) {
        // Clin-9/2012 determine rapidity more generally
        // C     to prevent overflow when Pt~=0 and E=|Pz|:
        // C            IF (ABS(PZAR(I)) .GE. PEAR(I)) THEN
        // C               PRINT *, ' IN ARINI1'
        // C               PRINT *, 'ABS(PZ) .GE. EE for particle ', I
        // C               PRINT *, ' FLAV = ', ITYPAR(I), ' PX = ', PXAR(I),
        // C     &              ' PY = ', PYAR(I)
        // C               PRINT *, ' PZ = ', PZAR(I), ' EE = ', PEAR(I)
        // C               PRINT *, ' XM = ', XMAR(I)
        // C               RAP = 1000000.0
        // C               GOTO 50
        // C            END IF
        // Cc            RAP=0.5*LOG((PEAR(I)+PZAR(I))/(PEAR(I)-PZAR(I)))
        // C RAP=0.5*LOG((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5)) C 50
        // CONTINUE
        if ((fem::pow2(xmar(i)) + fem::pow2(pxar(i)) + fem::pow2(pyar(i))) >
            0.f) {
          rap = asinh(pzar(i) /
                      fem::sqrt(fem::pow2(xmar(i)) + fem::pow2(pxar(i)) +
                                fem::pow2(pyar(i))));
        } else {
          write(6, star), " IN ARINI1 mt=0";
          rap = 1000000.0f * fem::sign(1.f, pzar(i));
        }
        // C
        vx = pxar(i) / pear(i);
        vy = pyar(i) / pear(i);
        ftar(i) = tau0 * fem::cosh(rap);
        gxar(i) += vx * ftar(i);
        gyar(i) += vy * ftar(i);
        gzar(i) = tau0 * fem::sinh(rap);
        // Clin-5/2009 No formation time for spectator projectile or target
        // nucleons:
        if (pxar(i) == 0 && pyar(i) == 0 &&
            (itypar(i) == 2112 || itypar(i) == 2212)) {
          // Clin-2/2013 for spectator target nucleons in LAB frame:
          // C     1           .and.(PEAR(I)*2/HINT1(1)).gt.0.99
          if ((pear(i) / hint1(6) > 0.99f && pear(i) / hint1(6) < 1.01f) ||
              (pear(i) / hint1(7) > 0.99f && pear(i) / hint1(7) < 1.01f)) {
            // C
            taui = 1.e-20f;
            ftar(i) = taui * fem::cosh(rap);
            gzar(i) = taui * fem::sinh(rap);
          }
        }
      }
      // Clin-8/2015:
    }
    // Clin-7/10/01-end
    // Clin-3/2009 cleanup of program flow:
  } else {
    FEM_DO_SAFE(i, 1, np) {
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZAR(I)) .GE. PEAR(I)) THEN
      // C               PRINT *, ' IN ARINI1'
      // C               PRINT *, 'ABS(PZ) .GE. EE for particle ', I
      // C               PRINT *, ' FLAV = ', ITYPAR(I), ' PX = ', PXAR(I),
      // C     &              ' PY = ', PYAR(I)
      // C               PRINT *, ' PZ = ', PZAR(I), ' EE = ', PEAR(I)
      // C               PRINT *, ' XM = ', XMAR(I)
      // C               RAP = 1000000.0
      // C               GOTO 100
      // Cc               STOP
      // C            END IF
      // C 100        CONTINUE
      // C            RAP=0.5*LOG((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
      if ((fem::pow2(xmar(i)) + fem::pow2(pxar(i)) + fem::pow2(pyar(i))) >
          0.f) {
        rap =
            asinh(pzar(i) / fem::sqrt(fem::pow2(xmar(i)) + fem::pow2(pxar(i)) +
                                      fem::pow2(pyar(i))));
      } else {
        write(6, star), " IN ARINI1 mt=0";
        rap = 1000000.0f * fem::sign(1.f, pzar(i));
      }
      // C
      vx = pxar(i) / pear(i);
      vy = pyar(i) / pear(i);
      // C.....give initial formation time shift
      taui = ftar(i) + tau0;
      ftar(i) = taui * fem::cosh(rap);
      gxar(i) += vx * tau0 * fem::cosh(rap);
      gyar(i) += vy * tau0 * fem::cosh(rap);
      // C     4/25/03: hadron z-position upon formation determined the same way
      // as x,y:
      gzar(i) = taui * fem::sinh(rap);
      // C     the old prescription:
      // C            GZAR(I) = GZAR(I) + TAU0 * SINH(RAP)
      zsmear = fem::sngl(cmn.smearh) * (2.f * ranart(cmn.nseed) - 1.f);
      gzar(i) += zsmear;
      // Cbz1/28/99end
      // C     10/05/01 no formation time for spectator projectile or target
      // nucleons:
      if (pxar(i) == 0 && pyar(i) == 0 &&
          (itypar(i) == 2112 || itypar(i) == 2212)) {
        // Clin-2/2013 for spectator target nucleons in LAB frame:
        // C     1           .and.(PEAR(I)*2/HINT1(1)).gt.0.99
        if ((pear(i) / hint1(6) > 0.99f && pear(i) / hint1(6) < 1.01f) ||
            (pear(i) / hint1(7) > 0.99f && pear(i) / hint1(7) < 1.01f)) {
          // C
          // Clin-5/2008:
          // C               TAUI=0.00001
          taui = 1.e-20f;
          ftar(i) = taui * fem::cosh(rap);
          gzar(i) = taui * fem::sinh(rap) + zsmear;
        }
      }
    }
    // Clin-3/2009 cleanup of program flow:
  }
  // C
  // Clin-3/2009 Add initial hadrons before the hadron cascade starts:
  addhad(cmn);
  // C
}

struct arindx_save {
  int i;
  int indxt;
  int ir;
  int j;
  int l;
  float q;

  arindx_save()
      : i(fem::int0),
        indxt(fem::int0),
        ir(fem::int0),
        j(fem::int0),
        l(fem::int0),
        q(fem::float0) {}
};

// C
// C=======================================================================
// C
// C.....Routine borrowed from ZPC.
// C.....double precision  is modified to real*4.
// C
// Cbz1/29/99
// C      subroutine index1(n, m, arrin, indx)
void arindx(common& cmn, int const& n, int const& m, arr_cref<float> arrin,
            arr_ref<int> indx) {
  FEM_CMN_SVE(arindx);
  arrin(dimension(n));
  indx(dimension(n));
  int& i = sve.i;
  int& indxt = sve.indxt;
  int& ir = sve.ir;
  int& j = sve.j;
  int& l = sve.l;
  float& q = sve.q;
  // Cbz1/29/99end
  // C     indexes the first m elements of ARRIN of length n, i.e., outputs INDX
  // C     such that ARRIN(INDEX(J)) is in ascending order for J=1,...,m
  // C
  // C      implicit real*4 (a-h, o-z)
  // C
  FEM_DO_SAFE(j, 1, m) { indx(j) = j; }
  l = m / 2 + 1;
  ir = m;
statement_10:
  if (l > 1) {
    l = l - 1;
    indxt = indx(l);
    q = arrin(indxt);
  } else {
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
    } else {
      j = ir + 1;
    }
    goto statement_20;
  }
  indx(i) = indxt;
  goto statement_10;
  // C
}

struct artord_save {
  static const int maxstr = 150001;

  arr<float> dptemp;
  arr<float> ee0;
  arr<float> ft0;
  arr<float> gx0;
  arr<float> gy0;
  arr<float> gz0;
  int i;
  arr<int> indx;
  arr<int> ityp0;
  int np;
  int npar;
  arr<float> px0;
  arr<float> py0;
  arr<float> pz0;
  arr<float> xm0;

  artord_save()
      : dptemp(dimension(maxstr), fem::fill0),
        ee0(dimension(maxstr), fem::fill0),
        ft0(dimension(maxstr), fem::fill0),
        gx0(dimension(maxstr), fem::fill0),
        gy0(dimension(maxstr), fem::fill0),
        gz0(dimension(maxstr), fem::fill0),
        i(fem::int0),
        indx(dimension(maxstr), fem::fill0),
        ityp0(dimension(maxstr), fem::fill0),
        np(fem::int0),
        npar(fem::int0),
        px0(dimension(maxstr), fem::fill0),
        py0(dimension(maxstr), fem::fill0),
        pz0(dimension(maxstr), fem::fill0),
        xm0(dimension(maxstr), fem::fill0) {}
};

const int artord_save::maxstr;

// C
// C-----------------------------------------------------------------------
// C
// C.....subroutine to order particle labels according to increasing
// C.....formation time
// C
void artord(common& cmn) {
  FEM_CMN_SVE(artord);
  // COMMON arprnt
  arr_ref<int> iaint2(cmn.iaint2, dimension(50));
  // COMMON arprc
  const int maxstr = 150001;
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
  // COMMON dpert
  arr_ref<float> dpertp(cmn.dpertp, dimension(maxstr));
  //
  // SAVE
  arr_ref<float> dptemp(sve.dptemp, dimension(maxstr));
  arr_ref<float> ee0(sve.ee0, dimension(maxstr));
  arr_ref<float> ft0(sve.ft0, dimension(maxstr));
  arr_ref<float> gx0(sve.gx0, dimension(maxstr));
  arr_ref<float> gy0(sve.gy0, dimension(maxstr));
  arr_ref<float> gz0(sve.gz0, dimension(maxstr));
  int& i = sve.i;
  arr_ref<int> indx(sve.indx, dimension(maxstr));
  arr_ref<int> ityp0(sve.ityp0, dimension(maxstr));
  int& np = sve.np;
  int& npar = sve.npar;
  arr_ref<float> px0(sve.px0, dimension(maxstr));
  arr_ref<float> py0(sve.py0, dimension(maxstr));
  arr_ref<float> pz0(sve.pz0, dimension(maxstr));
  arr_ref<float> xm0(sve.xm0, dimension(maxstr));
  //
  // C
  // C.....before invoking ARTORD:
  // C.....IAINT2(1) must be set:
  // Cc      SAVE /ARPRNT/
  // Cc      SAVE /ARPRC/
  // Clin-3/2009 Take care of particle weights when user inserts initial
  // hadrons:
  // C
  npar = 0;
  np = iaint2(1);
  FEM_DO_SAFE(i, 1, np) {
    ityp0(i) = itypar(i);
    gx0(i) = gxar(i);
    gy0(i) = gyar(i);
    gz0(i) = gzar(i);
    ft0(i) = ftar(i);
    px0(i) = pxar(i);
    py0(i) = pyar(i);
    pz0(i) = pzar(i);
    ee0(i) = pear(i);
    xm0(i) = xmar(i);
    // Clin-3/2009:
    dptemp(i) = dpertp(i);
  }
  arindx(cmn, maxstr, np, ft0, indx);
  FEM_DO_SAFE(i, 1, np) {
    // Cbz12/3/98
    // C         IF (ITYP0(INDX(I)) .EQ. 211) THEN
    // C         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 321) THEN
    // C         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 2212 .OR.
    // C     &      ITYP0(INDX(I)) .EQ. 2112 .OR. ITYP0(INDX(I)) .EQ. -211 .OR.
    // C     &      ITYP0(INDX(I)) .EQ. 111) THEN
    // C         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 2212 .OR.
    // C     &      ITYP0(INDX(I)) .EQ. 2112) THEN
    npar++;
    // C         ITYPAR(I) = ITYP0(INDX(I))
    // C         GXAR(I) = GX0(INDX(I))
    // C         GYAR(I) = GY0(INDX(I))
    // C         GZAR(I) = GZ0(INDX(I))
    // C         FTAR(I) = FT0(INDX(I))
    // C         PXAR(I) = PX0(INDX(I))
    // C         PYAR(I) = PY0(INDX(I))
    // C         PZAR(I) = PZ0(INDX(I))
    // C         PEAR(I) = EE0(INDX(I))
    // C         XMAR(I) = XM0(INDX(I))
    itypar(npar) = ityp0(indx(i));
    gxar(npar) = gx0(indx(i));
    gyar(npar) = gy0(indx(i));
    gzar(npar) = gz0(indx(i));
    ftar(npar) = ft0(indx(i));
    pxar(npar) = px0(indx(i));
    pyar(npar) = py0(indx(i));
    pzar(npar) = pz0(indx(i));
    pear(npar) = ee0(indx(i));
    xmar(npar) = xm0(indx(i));
    // Clin-3/2009:
    dpertp(npar) = dptemp(indx(i));
    // C         END IF
    // Cbz12/3/98end
  }
  iaint2(1) = npar;
  // C
}

struct arini_save {
  int iflg;

  arini_save() : iflg(fem::int0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....subroutine to initialize cascade.
// C
void arini(common& cmn) {
  FEM_CMN_SVE(arini);
  common_write write(cmn);
  arr_cref<int> iapar2(cmn.iapar2, dimension(50));
  //
  int& iflg = sve.iflg;
  // C
  // C.....before invoking ARINI:
  // C.....IAPAR2(1), IAINT2(1) must be set.
  // Cc      SAVE /ARPRNT/
  // C
  // Ctest off for resonance (phi, K*) studies:
  // C      OPEN (89, FILE = 'ana/decay_rec.dat', STATUS = 'UNKNOWN')
  // C
  iflg = iapar2(1);
  switch (iflg) {
    case 1:
      goto statement_200;
    case 2:
      goto statement_200;
    case 3:
      goto statement_300;
    default:
      break;
  }
  // C
  // C.....error choice of initialization
  write(6, star), "IAPAR2(1) must be 1, 2, or 3";
  FEM_STOP(0);
// C
// C.....to use default initial conditions generated by the cascade,
// C.....or to read in initial conditions.
statement_200:
  return;
// C
// C.....to generate formation time and the position at formation time from
// C.....read-in initial conditions with an averaged formation proper time.
statement_300:
  arini1(cmn);
  // C.....ordering the particle label according to increasing order of
  // C.....formation time.
  artord(cmn);
  // C
}

struct arini2_save {
  int i;
  int ip;
  int irun;

  arini2_save() : i(fem::int0), ip(fem::int0), irun(fem::int0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....subroutine to copy individually generated particle record into
// C.....particle record for many test particle runs.
// C
void arini2(common& cmn, int const& k) {
  FEM_CMN_SVE(arini2);
  // COMMON arprnt
  arr_cref<int> iaint2(cmn.iaint2, dimension(50));
  // COMMON arprc
  const int maxstr = 150001;
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
  // COMMON arerc1
  const int maxr = 1;
  arr_ref<int> multi1(cmn.multi1, dimension(maxr));
  // COMMON arprc1
  arr_ref<int, 2> ityp1(cmn.ityp1, dimension(maxstr, maxr));
  arr_ref<float, 2> gx1(cmn.gx1, dimension(maxstr, maxr));
  arr_ref<float, 2> gy1(cmn.gy1, dimension(maxstr, maxr));
  arr_ref<float, 2> gz1(cmn.gz1, dimension(maxstr, maxr));
  arr_ref<float, 2> ft1(cmn.ft1, dimension(maxstr, maxr));
  arr_ref<float, 2> px1(cmn.px1, dimension(maxstr, maxr));
  arr_ref<float, 2> py1(cmn.py1, dimension(maxstr, maxr));
  arr_ref<float, 2> pz1(cmn.pz1, dimension(maxstr, maxr));
  arr_ref<float, 2> ee1(cmn.ee1, dimension(maxstr, maxr));
  arr_ref<float, 2> xm1(cmn.xm1, dimension(maxstr, maxr));
  // COMMON tdecay
  arr_ref<float> tfdcy(cmn.tfdcy, dimension(maxstr));
  arr_ref<float, 2> tfdpi(cmn.tfdpi, dimension(maxstr, maxr));
  arr_ref<float> tft(cmn.tft, dimension(maxstr));
  // COMMON input1
  float& dt = cmn.dt;
  // COMMON input2
  int& ntmax = cmn.ntmax;
  // COMMON dpert
  arr_cref<float> dpertp(cmn.dpertp, dimension(maxstr));
  arr_ref<float, 2> dpp1(cmn.dpp1, dimension(maxstr, maxr));
  //
  // SAVE
  int& i = sve.i;
  int& ip = sve.ip;
  int& irun = sve.irun;
  //
  // C
  // Cc      SAVE /ARPRNT/
  // Cc      SAVE /ARPRC/
  // Cc      SAVE /ARERC1/
  // Cc      SAVE /ARPRC1/
  // Cc      SAVE /tdecay/
  // Cc      SAVE /input1/
  // Cc      SAVE /INPUT2/
  // Cc      SAVE /RNDF77/
  // C
  multi1(k) = iaint2(1);
  FEM_DO_SAFE(i, 1, multi1(k)) {
    ityp1(i, k) = itypar(i);
    gx1(i, k) = gxar(i);
    gy1(i, k) = gyar(i);
    gz1(i, k) = gzar(i);
    ft1(i, k) = ftar(i);
    px1(i, k) = pxar(i);
    py1(i, k) = pyar(i);
    pz1(i, k) = pzar(i);
    ee1(i, k) = pear(i);
    xm1(i, k) = xmar(i);
    // Clin-3/2009 hadron weights are initialized in addhad():
    // Clin-5/2008 all hadrons not perturbatively-produced have the weight of 1:
    // C         dpp1(I,K)=1.
    dpp1(i, k) = dpertp(i);
  }
  // C
  // C     initialize final time of each particle to ntmax*dt except for
  // C     decay daughters, which have values given by tfdcy() and >(ntmax*dt):
  FEM_DO_SAFE(ip, 1, maxstr) {
    tfdcy(ip) = ntmax * dt;
    tft(ip) = ntmax * dt;
  }
  // C
  FEM_DO_SAFE(irun, 1, maxr) {
    FEM_DO_SAFE(ip, 1, maxstr) { tfdpi(ip, irun) = ntmax * dt; }
  }
  // C
}

struct iarflv_save {
  float r;

  iarflv_save() : r(fem::float0) {}
};

// C
// C=======================================================================
// C
// C.....function to convert PDG flavor code into ART flavor code.
// C
int iarflv(common& cmn, int const& ipdg) {
  int return_value = fem::int0;
  FEM_CMN_SVE(iarflv);
  // SAVE
  float& r = sve.r;
  //
  // C
  // Cc      SAVE /input1/
  // Cc      SAVE /RNDF77/
  // C
  // C.....anti-Delta-
  if (ipdg == -1114) {
    return_value = -6;
    return return_value;
  }
  // C
  // C.....anti-Delta0
  if (ipdg == -2114) {
    return_value = -7;
    return return_value;
  }
  // C
  // C.....anti-Delta+
  if (ipdg == -2214) {
    return_value = -8;
    return return_value;
  }
  // C
  // C.....anti-Delta++
  if (ipdg == -2224) {
    return_value = -9;
    return return_value;
  }
  // C
  // Cbzdbg2/23/99
  // C.....anti-proton
  if (ipdg == -2212) {
    return_value = -1;
    return return_value;
  }
  // C
  // C.....anti-neutron
  if (ipdg == -2112) {
    return_value = -2;
    return return_value;
  }
  // Cbzdbg2/23/99end
  // C
  // C.....eta
  if (ipdg == 221) {
    return_value = 0;
    return return_value;
  }
  // C
  // C.....proton
  if (ipdg == 2212) {
    return_value = 1;
    return return_value;
  }
  // C
  // C.....neutron
  if (ipdg == 2112) {
    return_value = 2;
    return return_value;
  }
  // C
  // C.....pi-
  if (ipdg == -211) {
    return_value = 3;
    return return_value;
  }
  // C
  // C.....pi0
  if (ipdg == 111) {
    return_value = 4;
    return return_value;
  }
  // C
  // C.....pi+
  if (ipdg == 211) {
    return_value = 5;
    return return_value;
  }
  // C
  // C.....Delta-
  if (ipdg == 1114) {
    return_value = 6;
    return return_value;
  }
  // C
  // C.....Delta0
  if (ipdg == 2114) {
    return_value = 7;
    return return_value;
  }
  // C
  // C.....Delta+
  if (ipdg == 2214) {
    return_value = 8;
    return return_value;
  }
  // C
  // C.....Delta++
  if (ipdg == 2224) {
    return_value = 9;
    return return_value;
  }
  // C
  // C.....Lambda
  if (ipdg == 3122) {
    return_value = 14;
    return return_value;
  }
  // C
  // C.....Lambda-bar
  if (ipdg == -3122) {
    return_value = -14;
    return return_value;
  }
  // C
  // C.....Sigma-
  if (ipdg == 3112) {
    return_value = 15;
    return return_value;
  }
  // C
  // C.....Sigma-bar
  if (ipdg == -3112) {
    return_value = -15;
    return return_value;
  }
  // C
  // C.....Sigma0
  if (ipdg == 3212) {
    return_value = 16;
    return return_value;
  }
  // C
  // C.....Sigma0-bar
  if (ipdg == -3212) {
    return_value = -16;
    return return_value;
  }
  // C
  // C.....Sigma+
  if (ipdg == 3222) {
    return_value = 17;
    return return_value;
  }
  // C
  // C.....Sigma+ -bar
  if (ipdg == -3222) {
    return_value = -17;
    return return_value;
  }
  // C
  // C.....K-
  if (ipdg == -321) {
    return_value = 21;
    return return_value;
  }
  // C
  // C.....K+
  if (ipdg == 321) {
    return_value = 23;
    return return_value;
  }
  // C
  // C.....temporary entry for K0
  if (ipdg == 311) {
    return_value = 23;
    return return_value;
  }
  // C
  // C.....temporary entry for K0bar
  if (ipdg == -311) {
    return_value = 21;
    return return_value;
  }
  // C
  // C.....temporary entry for K0S and K0L
  if (ipdg == 310 || ipdg == 130) {
    r = ranart(cmn.nseed);
    if (r > 0.5f) {
      return_value = 23;
    } else {
      return_value = 21;
    }
    return return_value;
  }
  // C
  // C.....rho-
  if (ipdg == -213) {
    return_value = 25;
    return return_value;
  }
  // C
  // C.....rho0
  if (ipdg == 113) {
    return_value = 26;
    return return_value;
  }
  // C
  // C.....rho+
  if (ipdg == 213) {
    return_value = 27;
    return return_value;
  }
  // C
  // C.....omega
  if (ipdg == 223) {
    return_value = 28;
    return return_value;
  }
  // C
  // C.....phi
  if (ipdg == 333) {
    return_value = 29;
    return return_value;
  }
  // C
  // C.....K*+
  if (ipdg == 323) {
    return_value = 30;
    return return_value;
  }
  // C.....K*-
  if (ipdg == -323) {
    return_value = -30;
    return return_value;
  }
  // C.....temporary entry for K*0
  if (ipdg == 313) {
    return_value = 30;
    return return_value;
  }
  // C.....temporary entry for K*0bar
  if (ipdg == -313) {
    return_value = -30;
    return return_value;
  }
  // C
  // C...... eta-prime
  if (ipdg == 331) {
    return_value = 31;
    return return_value;
  }
  // C
  // C...... a1
  // C     IF (IPDG .EQ. 777) THEN
  // C        IARFLV = 32
  // C        RETURN
  // C     END IF
  // C
  // C... cascade-
  if (ipdg == 3312) {
    return_value = 40;
    return return_value;
  }
  // C
  // C... cascade+ (bar)
  if (ipdg == -3312) {
    return_value = -40;
    return return_value;
  }
  // C
  // C... cascade0
  if (ipdg == 3322) {
    return_value = 41;
    return return_value;
  }
  // C
  // C... cascade0 -bar
  if (ipdg == -3322) {
    return_value = -41;
    return return_value;
  }
  // C
  // C... Omega-
  if (ipdg == 3334) {
    return_value = 45;
    return return_value;
  }
  // C
  // C... Omega+ (bar)
  if (ipdg == -3334) {
    return_value = -45;
    return return_value;
  }
  // C
  // C... Di-Omega
  if (ipdg == 6666) {
    return_value = 44;
    return return_value;
  }
  // C sp06/05/01 end
  // C
  // Clin-3/2009 keep the same ID numbers in case there are initial deuterons:
  if (ipdg == 42 || ipdg == -42) {
    return_value = ipdg;
    return return_value;
  }
  // C
  // C.....other
  return_value = ipdg + 10000;
  // C
  return return_value;
}

// C
// C-----------------------------------------------------------------------
// C
// C.....function to convert ART flavor code into PDG flavor code.
// C
int invflv(common& cmn, int const& iart) {
  int return_value = fem::int0;
  // COMMON rndf77
  int& nseed = cmn.nseed;
  //
  // C
  // Cc      SAVE /input1/
  // Cc      SAVE /RNDF77/
  // C
  // C.....anti-Delta-
  if (iart == -6) {
    return_value = -1114;
    return return_value;
  }
  // C
  // C.....anti-Delta0
  if (iart == -7) {
    return_value = -2114;
    return return_value;
  }
  // C
  // C.....anti-Delta+
  if (iart == -8) {
    return_value = -2214;
    return return_value;
  }
  // C
  // C.....anti-Delta++
  if (iart == -9) {
    return_value = -2224;
    return return_value;
  }
  // C
  // Cbzdbg2/23/99
  // C.....anti-proton
  if (iart == -1) {
    return_value = -2212;
    return return_value;
  }
  // C
  // C.....anti-neutron
  if (iart == -2) {
    return_value = -2112;
    return return_value;
  }
  // Cbzdbg2/23/99end
  // C
  // C.....eta
  if (iart == 0) {
    return_value = 221;
    return return_value;
  }
  // C
  // C.....proton
  if (iart == 1) {
    return_value = 2212;
    return return_value;
  }
  // C
  // C.....neutron
  if (iart == 2) {
    return_value = 2112;
    return return_value;
  }
  // C
  // C.....pi-
  if (iart == 3) {
    return_value = -211;
    return return_value;
  }
  // C
  // C.....pi0
  if (iart == 4) {
    return_value = 111;
    return return_value;
  }
  // C
  // C.....pi+
  if (iart == 5) {
    return_value = 211;
    return return_value;
  }
  // C
  // C.....Delta-
  if (iart == 6) {
    return_value = 1114;
    return return_value;
  }
  // C
  // C.....Delta0
  if (iart == 7) {
    return_value = 2114;
    return return_value;
  }
  // C
  // C.....Delta+
  if (iart == 8) {
    return_value = 2214;
    return return_value;
  }
  // C
  // C.....Delta++
  if (iart == 9) {
    return_value = 2224;
    return return_value;
  }
  // C
  // Cc.....N*(1440), N*(1535) temporary entry
  // C      IF (IART .GE. 10 .AND. IART .LE.13) THEN
  // C         INVFLV = 0
  // C         RETURN
  // C      END IF
  // C
  // C.....Lambda
  if (iart == 14) {
    return_value = 3122;
    return return_value;
  }
  // C.....Lambda-bar
  if (iart == -14) {
    return_value = -3122;
    return return_value;
  }
  // C
  // Cbz3/12/99
  // C.....temporary entry for Sigma's
  // C      IF (IART .EQ. 15) THEN
  // C         R = RANART(NSEED)
  // C         IF (R .GT. 2. / 3.) THEN
  // C            INVFLV = 3112
  // C         ELSE IF (R .GT. 1./ 3. .AND. R .LE. 2. / 3.) THEN
  // C            INVFLV = 3212
  // C         ELSE
  // C            INVFLV = 3222
  // C         END IF
  // C         RETURN
  // C      END IF
  // C
  // C.....Sigma-
  if (iart == 15) {
    return_value = 3112;
    return return_value;
  }
  // C
  // C.....Sigma- bar
  if (iart == -15) {
    return_value = -3112;
    return return_value;
  }
  // C
  // C.....Sigma0
  if (iart == 16) {
    return_value = 3212;
    return return_value;
  }
  // C
  // C.....Sigma0 -bar
  if (iart == -16) {
    return_value = -3212;
    return return_value;
  }
  // C
  // C.....Sigma+
  if (iart == 17) {
    return_value = 3222;
    return return_value;
  }
  // C
  // C.....Sigma+ -bar
  if (iart == -17) {
    return_value = -3222;
    return return_value;
  }
  // C
  // Clin-2/23/03 K0S and K0L are generated at the last timestep:
  // C.....temporary entry for K- and K0bar
  if (iart == 21) {
    // C         R = RANART(NSEED)
    // C         IF (R .GT. 0.5) THEN
    return_value = -321;
    // C         ELSE
    // C            INVFLV = -311
    // C            R = RANART(NSEED)
    // C            IF (R .GT. 0.5) THEN
    // C               INVFLV = 310
    // C            ELSE
    // C               INVFLV = 130
    // C            END IF
    // C         END IF
    return return_value;
  }
  // C
  // C.....temporary entry for K+ and K0
  if (iart == 23) {
    // C         R = RANART(NSEED)
    // C         IF (R .GT. 0.5) THEN
    return_value = 321;
    // C         ELSE
    // C            INVFLV = 311
    // C            R = RANART(NSEED)
    // C            IF (R .GT. 0.5) THEN
    // C               INVFLV = 310
    // C            ELSE
    // C               INVFLV = 130
    // C            END IF
    // C         END IF
    return return_value;
  }
  // C
  // C.....K0Long:
  if (iart == 22) {
    return_value = 130;
    return return_value;
  }
  // C.....K0Short:
  if (iart == 24) {
    return_value = 310;
    return return_value;
  }
  // C
  // C.....rho-
  if (iart == 25) {
    return_value = -213;
    return return_value;
  }
  // C
  // C.....rho0
  if (iart == 26) {
    return_value = 113;
    return return_value;
  }
  // C
  // C.....rho+
  if (iart == 27) {
    return_value = 213;
    return return_value;
  }
  // C
  // C.....omega
  if (iart == 28) {
    return_value = 223;
    return return_value;
  }
  // C
  // C.....phi
  if (iart == 29) {
    return_value = 333;
    return return_value;
  }
  // C
  // C.....temporary entry for K*+ and K*0
  if (iart == 30) {
    return_value = 323;
    if (ranart(nseed) > 0.5f) {
      return_value = 313;
    }
    return return_value;
  }
  // C
  // C.....temporary entry for K*- and K*0bar
  if (iart == -30) {
    return_value = -323;
    if (ranart(nseed) > 0.5f) {
      return_value = -313;
    }
    return return_value;
  }
  // C
  // C... eta-prime (bar)
  if (iart == 31) {
    return_value = 331;
    return return_value;
  }
  // C
  // C... a1
  if (iart == 32) {
    return_value = 777;
    return return_value;
  }
  // C
  // C... cascade-
  if (iart == 40) {
    return_value = 3312;
    return return_value;
  }
  // C
  // C... cascade+ (bar)
  if (iart == -40) {
    return_value = -3312;
    return return_value;
  }
  // C
  // C... cascade0
  if (iart == 41) {
    return_value = 3322;
    return return_value;
  }
  // C
  // C... cascade0 -bar
  if (iart == -41) {
    return_value = -3322;
    return return_value;
  }
  // C
  // C... Omega-
  if (iart == 45) {
    return_value = 3334;
    return return_value;
  }
  // C
  // C... Omega+ (bar)
  if (iart == -45) {
    return_value = -3334;
    return return_value;
  }
  // C
  // C... Di-Omega
  if (iart == 44) {
    return_value = 6666;
    return return_value;
  }
  // C sp 12/19/00 end
  // C
  // Clin-5/2008 deuteron ID numbers in ART and ampt.dat:
  if (iart == 42) {
    return_value = 42;
    return return_value;
  } else if (iart == -42) {
    return_value = -42;
    return return_value;
  }
  // C
  // C.....other
  return_value = iart - 10000;
  // C
  return return_value;
}

struct blockdata_ardata_save {};

// C
// C=======================================================================
// C
void blockdata_ardata(common& cmn) {
  FEM_CMN_SVE(blockdata_ardata);
  // COMMON arprnt
  arr_ref<float> arpar1(cmn.arpar1, dimension(100));
  arr_ref<int> iapar2(cmn.iapar2, dimension(50));
  arr_ref<float> arint1(cmn.arint1, dimension(100));
  arr_ref<int> iaint2(cmn.iaint2, dimension(50));
  //
  if (is_called_first_time) {
    fem::data((values, 1.19f, 99 * datum(0.0f))), arpar1;
    fem::data((values, 3, 49 * datum(0))), iapar2;
    fem::data((values, 100 * datum(0.0f))), arint1;
    fem::data((values, 50 * datum(0))), iaint2;
  }
  // C
  // Cc      SAVE /ARPRNT/
  // C
}

struct x2kaon_save {
  float f1;
  float f2;
  float f3;
  float sigma1;
  float sigma2;
  float sigma3;
  float smin;
  float x;

  x2kaon_save()
      : f1(fem::float0),
        f2(fem::float0),
        f3(fem::float0),
        sigma1(fem::float0),
        sigma2(fem::float0),
        sigma3(fem::float0),
        smin(fem::float0),
        x(fem::float0) {}
};

// C
// C*****************************************
// C for pp-->pp + kaon + anti-kaon
// C      real*4 function X2kaon(srt)
float x2kaon(common& cmn, float const& srt) {
  float return_value = fem::float0;
  FEM_CMN_SVE(x2kaon);
  // SAVE
  float& f1 = sve.f1;
  float& f2 = sve.f2;
  float& f3 = sve.f3;
  float& sigma1 = sve.sigma1;
  float& sigma2 = sve.sigma2;
  float& sigma3 = sve.sigma3;
  float& smin = sve.smin;
  float& x = sve.x;
  //
  // C  This function contains the experimental total pp->pp+K(+)K(-) Xsections
  // * C  srt    = DSQRT(s) in GeV * C  xsec   = production cross section in mb *
  // C * C***************************************** C     minimum c.m.s. energy
  // to create 2 kaon. = 2*(mp+mk)
  smin = 2.8639f;
  return_value = 0.0000001f;
  if (srt < smin) {
    return return_value;
  }
  sigma1 = 2.8f;
  sigma2 = 7.7f;
  sigma3 = 3.9f;
  x = fem::pow2(srt) / fem::pow2(smin) + 0.0000001f;
  f1 = (1.f + 1.f / fem::sqrt(x)) * fem::alog(x) -
       4.f * (1.f - 1.f / fem::sqrt(x));
  f2 = 1.f - (1.f / fem::sqrt(x)) * (1.f + fem::alog(fem::sqrt(x)));
  f3 = fem::pow(((x - 1.f) / fem::pow2(x)), 3.5f);
  return_value =
      fem::pow3((1.f - 1.f / x)) * (sigma1 * f1 + sigma2 * f2) + sigma3 * f3;
  return return_value;
}

struct pinsg0_save {
  float ratio;
  float srt0;

  pinsg0_save() : ratio(fem::float0), srt0(fem::float0) {}
};

float pinsg0(common& cmn, float const& srt) {
  float return_value = fem::float0;
  FEM_CMN_SVE(pinsg0);
  // SAVE
  float& ratio = sve.ratio;
  float& srt0 = sve.srt0;
  //
  // C cross section in mb for PI- + P -> P + K0 + K-
  // C     Mn + 2* Mk
  srt0 = 0.938f + 2.f * 0.498f;
  if (srt < srt0) {
    return_value = 0.0f;
    return return_value;
  }
  ratio = fem::pow2(srt0) / fem::pow2(srt);
  return_value = 1.121f * fem::pow((1.f - ratio), 1.86f) * fem::pow2(ratio);
  return return_value;
}

struct aknel_save {
  float sigma1;

  aknel_save() : sigma1(fem::float0) {}
};

float aknel(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  FEM_CMN_SVE(aknel);
  // SAVE
  float& sigma1 = sve.sigma1;
  //
  // Ccross section in mb for K- + N reactions.
  // C        the following data come from PRC 41 (1701)
  // C        sigma1: K(-) + neutron elastic
  if (pkaon < 0.5f || pkaon >= 4.0f) {
    sigma1 = 0.f;
  }
  if (pkaon >= 0.5f && pkaon < 1.0f) {
    sigma1 = 20.f * fem::pow(pkaon, 2.74f);
  }
  if (pkaon >= 1.0f && pkaon < 4.0f) {
    sigma1 = 20.f * fem::pow(pkaon, (-1.8f));
  }
  return_value = sigma1;
  return return_value;
}

struct akpel_save {
  float sigma2;

  akpel_save() : sigma2(fem::float0) {}
};

float akpel(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  FEM_CMN_SVE(akpel);
  // SAVE
  float& sigma2 = sve.sigma2;
  //
  // Ccross section in mb for K- + N reactions.
  // C        the following data come from PRC 41 (1701)
  // C        sigma2: K(-) + proton elastic
  if (pkaon < 0.25f || pkaon >= 4.0f) {
    sigma2 = 0.f;
  }
  if (pkaon >= 0.25f && pkaon < 4.0f) {
    sigma2 = 13.f * fem::pow(pkaon, (-0.9f));
  }
  return_value = sigma2;
  return return_value;
}

struct aknsgm_save {
  float sigma2;

  aknsgm_save() : sigma2(fem::float0) {}
};

float aknsgm(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  FEM_CMN_SVE(aknsgm);
  // SAVE
  float& sigma2 = sve.sigma2;
  //
  // Ccross section in mb for K- + N reactions.
  // C        sigma2: x section for K- + n -> sigma0 + PI-
  if (pkaon < 0.5f || pkaon >= 6.0f) {
    sigma2 = 0.f;
  }
  if (pkaon >= 0.5f && pkaon < 1.0f) {
    sigma2 = 1.2f * fem::pow(pkaon, (-1.3f));
  }
  if (pkaon >= 1.0f && pkaon < 6.0f) {
    sigma2 = 1.2f * fem::pow(pkaon, (-2.3f));
  }
  return_value = sigma2;
  return return_value;
}

struct akpsgm_save {
  float sigma1;

  akpsgm_save() : sigma1(fem::float0) {}
};

float akpsgm(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  FEM_CMN_SVE(akpsgm);
  // SAVE
  float& sigma1 = sve.sigma1;
  //
  // Ccross section in mb for K- + N reactions.
  // C        sigma1: x section for K- + p -> sigma0 + PI0
  if (pkaon < 0.2f || pkaon >= 1.5f) {
    sigma1 = 0.f;
  }
  if (pkaon >= 0.2f && pkaon < 1.5f) {
    sigma1 = 0.6f * fem::pow(pkaon, (-1.8f));
  }
  return_value = sigma1;
  return return_value;
}

struct akplam_save {
  float p;
  float sigma;

  akplam_save() : p(fem::float0), sigma(fem::float0) {}
};

float akplam(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  FEM_CMN_SVE(akplam);
  // SAVE
  float& p = sve.p;
  float& sigma = sve.sigma;
  //
  // Ccross section in mb for K- + N reactions.
  // C        sigma: x section for K- + p -> lambda + PI0
  p = pkaon;
  if (pkaon < 0.2f || pkaon >= 10.0f) {
    sigma = 0.f;
  }
  if (pkaon >= 0.2f && pkaon < 0.9f) {
    sigma = 50.f * fem::pow2(p) - 67.f * p + 24.f;
  }
  if (pkaon >= 0.9f && pkaon < 10.0f) {
    sigma = 3.0f * fem::pow(pkaon, (-2.6f));
  }
  return_value = sigma;
  return return_value;
}

float aknlam(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  // Ccross section in mb for K- + N reactions.
  return_value = akplam(cmn, pkaon);
  return return_value;
}

struct fstate_save {
  float aka;
  float bbb;
  float beta;
  float ekm;
  float ekmax;
  float ekp;
  float fac;
  float gama;
  float guass;
  int icount;
  arr<float> pe;
  float pio;
  float pkmax;
  float ptkm;
  float ptkmi2;
  float ptkp;
  float ptkpl2;
  float ptp;
  float ptp2;
  float pzcms;
  float pzkm;
  float pzkp;
  float resten;
  float restms;
  float restpz;
  float rmt3;
  float rmt4;
  float rsq;
  float v1;
  float v2;
  float xstar;

  fstate_save()
      : aka(fem::float0),
        bbb(fem::float0),
        beta(fem::float0),
        ekm(fem::float0),
        ekmax(fem::float0),
        ekp(fem::float0),
        fac(fem::float0),
        gama(fem::float0),
        guass(fem::float0),
        icount(fem::int0),
        pe(dimension(4), fem::fill0),
        pio(fem::float0),
        pkmax(fem::float0),
        ptkm(fem::float0),
        ptkmi2(fem::float0),
        ptkp(fem::float0),
        ptkpl2(fem::float0),
        ptp(fem::float0),
        ptp2(fem::float0),
        pzcms(fem::float0),
        pzkm(fem::float0),
        pzkp(fem::float0),
        resten(fem::float0),
        restms(fem::float0),
        restpz(fem::float0),
        rmt3(fem::float0),
        rmt4(fem::float0),
        rsq(fem::float0),
        v1(fem::float0),
        v2(fem::float0),
        xstar(fem::float0) {}
};

void fstate(common& cmn, int const& /* iseed */, float const& srt,
            float const& dm3, float const& dm4, arr_ref<float> px,
            arr_ref<float> py, arr_ref<float> pz, int& iflag) {
  FEM_CMN_SVE(fstate);
  px(dimension(4));
  py(dimension(4));
  pz(dimension(4));
  common_write write(cmn);
  int& nseed = cmn.nseed;
  //
  float& aka = sve.aka;
  float& bbb = sve.bbb;
  float& beta = sve.beta;
  float& ekm = sve.ekm;
  float& ekmax = sve.ekmax;
  float& ekp = sve.ekp;
  float& fac = sve.fac;
  float& gama = sve.gama;
  float& guass = sve.guass;
  int& icount = sve.icount;
  arr_ref<float> pe(sve.pe, dimension(4));
  float& pio = sve.pio;
  float& pkmax = sve.pkmax;
  float& ptkm = sve.ptkm;
  float& ptkmi2 = sve.ptkmi2;
  float& ptkp = sve.ptkp;
  float& ptkpl2 = sve.ptkpl2;
  float& ptp = sve.ptp;
  float& ptp2 = sve.ptp2;
  float& pzcms = sve.pzcms;
  float& pzkm = sve.pzkm;
  float& pzkp = sve.pzkp;
  float& resten = sve.resten;
  float& restms = sve.restms;
  float& restpz = sve.restpz;
  float& rmt3 = sve.rmt3;
  float& rmt4 = sve.rmt4;
  float& rsq = sve.rsq;
  float& v1 = sve.v1;
  float& v2 = sve.v2;
  float& xstar = sve.xstar;
  // C        function: decide final momentum for N,N,K(+),and K(-)
  // Cc      SAVE /RNDF77/
  // C
  iflag = -1;
  // C        iflag=-1: fail to find momenta
  // C             = 1: success
  pio = 3.1415926f;
  aka = 0.498f;
  // C        v=0.43
  // C        w=-0.84
  // C        b=3.78
  // C        c=0.47
  // C        d=3.60
  // C        fmax=1.056
  // C        gmax=1.+c
  // C
  icount = 0;
  ekmax = (srt - dm3 - dm4) / 2.f;
  if (ekmax <= aka) {
    return;
  }
  pkmax = fem::sqrt(fem::pow2(ekmax) - fem::pow2(aka));
  // C
  if (dm3 <= 0.0f || dm4 <= 0.0f) {
    write(1, star), "error: minus mass!!!";
    return;
  }
// C
// C        after we have the momenta for both nucleus, we sample the
// C        transverse momentum for K-.
// C        dsigma/dpt**2 = exp(-4.145*pt**2) obtained by fitting data on
// C        page 72, fig 23i.
statement_50:
  icount++;
  if (icount > 10) {
    return;
  }
  ptkmi2 = -1.f / 4.145f * fem::alog(ranart(nseed));
  ptkm = fem::sqrt(ptkmi2);
statement_3:
  v1 = ranart(nseed);
  v2 = ranart(nseed);
  rsq = fem::pow2(v1) + fem::pow2(v2);
  if (rsq >= 1.0f || rsq <= 0.f) {
    goto statement_3;
  }
  fac = fem::sqrt(-2.f * fem::alog(rsq) / rsq);
  guass = v1 * fac;
  if (guass >= 5.f) {
    goto statement_3;
  }
  xstar = guass / 5.f;
  pzkm = pkmax * xstar;
  ekm = fem::sqrt(fem::pow2(aka) + fem::pow2(pzkm) + fem::pow2(ptkm));
  if (ranart(nseed) > aka / ekm) {
    goto statement_50;
  }
  bbb = ranart(nseed);
  px(3) = ptkm * fem::cos(2.f * pio * bbb);
  py(3) = ptkm * fem::sin(2.f * pio * bbb);
  if (ranart(nseed) > 0.5f) {
    pzkm = -1.f * pzkm;
  }
  pz(3) = pzkm;
  pe(3) = ekm;
statement_150:
  ptkpl2 = -1.f / 3.68f * fem::alog(ranart(nseed));
  ptkp = fem::sqrt(ptkpl2);
statement_13:
  v1 = ranart(nseed);
  v2 = ranart(nseed);
  rsq = fem::pow2(v1) + fem::pow2(v2);
  if (rsq >= 1.0f || rsq <= 0.f) {
    goto statement_13;
  }
  fac = fem::sqrt(-2.f * fem::alog(rsq) / rsq);
  guass = v1 * fac;
  if (guass >= 3.25f) {
    goto statement_13;
  }
  xstar = guass / 3.25f;
  pzkp = pkmax * xstar;
  ekp = fem::sqrt(fem::pow2(aka) + fem::pow2(pzkp) + fem::pow2(ptkp));
  if (ranart(nseed) > aka / ekp) {
    goto statement_150;
  }
  bbb = ranart(nseed);
  px(4) = ptkp * fem::cos(2.f * pio * bbb);
  py(4) = ptkp * fem::sin(2.f * pio * bbb);
  if (ranart(nseed) > 0.5f) {
    pzkp = -1.f * pzkp;
  }
  pz(4) = pzkp;
  pe(4) = ekp;
  // C
  resten = srt - pe(3) - pe(4);
  restpz = -pz(3) - pz(4);
  // C     resample
  if (resten <= fem::abs(restpz)) {
    goto statement_50;
  }
  restms = fem::sqrt(fem::pow2(resten) - fem::pow2(restpz));
  // C     resample
  if (restms < (dm3 + dm4)) {
    goto statement_50;
  }
  ptp2 = -1.f / 2.76f * fem::alog(ranart(nseed));
  ptp = fem::sqrt(ptp2);
  bbb = ranart(nseed);
  px(2) = ptp * fem::cos(2.f * pio * bbb);
  py(2) = ptp * fem::sin(2.f * pio * bbb);
  px(1) = -1.f * (px(4) + px(3) + px(2));
  py(1) = -1.f * (py(4) + py(3) + py(2));
  // C     transverse mass for K-
  rmt3 = fem::sqrt(fem::pow2(dm3) + fem::pow2(px(1)) + fem::pow2(py(1)));
  // C     transverse mass for K+
  rmt4 = fem::sqrt(fem::pow2(dm4) + fem::pow2(px(2)) + fem::pow2(py(2)));
  if (restms < (rmt3 + rmt4)) {
    goto statement_50;
  }
  // C        else: sampling success!
  pzcms = fem::sqrt((fem::pow2(restms) - fem::pow2((rmt3 + rmt4))) *
                    (fem::pow2(restms) - fem::pow2((rmt3 - rmt4)))) /
          2.f / restms;
  if (ranart(nseed) > 0.5f) {
    pz(1) = pzcms;
    pz(2) = -pzcms;
  } else {
    pz(1) = -pzcms;
    pz(2) = pzcms;
  }
  beta = restpz / resten;
  gama = 1.f / fem::sqrt(1.f - fem::pow2(beta));
  pz(1) = pz(1) * gama +
          beta * gama * fem::sqrt(fem::pow2(rmt3) + fem::pow2(pz(1)));
  pz(2) = pz(2) * gama +
          beta * gama * fem::sqrt(fem::pow2(rmt4) + fem::pow2(pz(2)));
  pe(1) = fem::sqrt(fem::pow2(rmt3) + fem::pow2(pz(1)));
  pe(2) = fem::sqrt(fem::pow2(rmt4) + fem::pow2(pz(2)));
  // C
  iflag = 1;
}

struct nnkaon_save {
  float betaak;
  float betak;
  float dm3;
  float dm4;
  float e1cm;
  float e2cm;
  float epcmak;
  float epcmk;
  float eti1;
  float eti2;
  int iflag;
  int lb1;
  int lb2;
  int n;
  float p1beta;
  float p2beta;
  float pt1i1;
  float pt1i2;
  float pt2i1;
  float pt2i2;
  float pt3i1;
  float pt3i2;
  arr<float> px;
  float pxrota;
  arr<float> py;
  float pyrota;
  arr<float> pz;
  float pzrota;
  float transf;

  nnkaon_save()
      : betaak(fem::float0),
        betak(fem::float0),
        dm3(fem::float0),
        dm4(fem::float0),
        e1cm(fem::float0),
        e2cm(fem::float0),
        epcmak(fem::float0),
        epcmk(fem::float0),
        eti1(fem::float0),
        eti2(fem::float0),
        iflag(fem::int0),
        lb1(fem::int0),
        lb2(fem::int0),
        n(fem::int0),
        p1beta(fem::float0),
        p2beta(fem::float0),
        pt1i1(fem::float0),
        pt1i2(fem::float0),
        pt2i1(fem::float0),
        pt2i2(fem::float0),
        pt3i1(fem::float0),
        pt3i2(fem::float0),
        px(dimension(4), fem::fill0),
        pxrota(fem::float0),
        py(dimension(4), fem::fill0),
        pyrota(fem::float0),
        pz(dimension(4), fem::fill0),
        pzrota(fem::float0),
        transf(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....extracted from G. Song's ART expasion including K- interactions
// C.....file `NEWNNK.FOR'
// C
void nnkaon(common& cmn, int const& irun, int const& iseed, int& ictrl,
            int const& i1, int const& i2, int& iblock, float const& srt,
            float const& pcx, float const& pcy, float const& pcz,
            int const& nchrg) {
  FEM_CMN_SVE(nnkaon);
  // COMMON aa
  const int maxstr = 150001;
  arr_cref<float, 2> r(cmn.r, dimension(3, maxstr));
  // COMMON bb
  arr_ref<float, 2> p(cmn.p, dimension(3, maxstr));
  // COMMON cc
  arr_ref<float> e(cmn.e, dimension(maxstr));
  // COMMON ee
  arr_ref<int> lb(cmn.lb, dimension(maxstr));
  // COMMON bg
  float& betax = cmn.betax;
  float& betay = cmn.betay;
  float& betaz = cmn.betaz;
  float& gamma = cmn.gamma;
  // COMMON nn
  int& nnn = cmn.nnn;
  // COMMON pa
  const int maxr = 1;
  arr_ref<float, 3> rpion(cmn.rpion, dimension(3, maxstr, maxr));
  // COMMON pb
  arr_ref<float, 3> ppion(cmn.ppion, dimension(3, maxstr, maxr));
  // COMMON pc
  arr_ref<float, 2> epion(cmn.epion, dimension(maxstr, maxr));
  // COMMON pd
  arr_ref<int, 2> lpion(cmn.lpion, dimension(maxstr, maxr));
  // COMMON dpert
  arr_cref<float> dpertp(cmn.dpertp, dimension(maxstr));
  arr_ref<float, 2> dppion(cmn.dppion, dimension(maxstr, maxr));
  //
  // SAVE
  float& betaak = sve.betaak;
  float& betak = sve.betak;
  float& dm3 = sve.dm3;
  float& dm4 = sve.dm4;
  float& e1cm = sve.e1cm;
  float& e2cm = sve.e2cm;
  float& epcmak = sve.epcmak;
  float& epcmk = sve.epcmk;
  float& eti1 = sve.eti1;
  float& eti2 = sve.eti2;
  int& iflag = sve.iflag;
  int& lb1 = sve.lb1;
  int& lb2 = sve.lb2;
  int& n = sve.n;
  float& p1beta = sve.p1beta;
  float& p2beta = sve.p2beta;
  float& pt1i1 = sve.pt1i1;
  float& pt1i2 = sve.pt1i2;
  float& pt2i1 = sve.pt2i1;
  float& pt2i2 = sve.pt2i2;
  float& pt3i1 = sve.pt3i1;
  float& pt3i2 = sve.pt3i2;
  arr_ref<float> px(sve.px, dimension(4));
  float& pxrota = sve.pxrota;
  arr_ref<float> py(sve.py, dimension(4));
  float& pyrota = sve.pyrota;
  arr_ref<float> pz(sve.pz, dimension(4));
  float& pzrota = sve.pzrota;
  float& transf = sve.transf;
  //
  // C        <pt>=0.27+0.037*log(srt) was changed to 0.632 + ... on Aug. 14,
  // 1997 C     CANCELED also alpha=1 changed to alpha=3 to decrease the leadng
  // effect. Cc      SAVE /AA/ Cc      SAVE /BB/ Cc      SAVE /CC/ Cc      SAVE
  // /EE/ Cc      SAVE /BG/ Cc      SAVE /NN/ Cc      SAVE /RUN/ Cc      SAVE
  // /PA/ Cc      SAVE /PB/ Cc      SAVE /PC/ Cc      SAVE /PD/ C      dm1=e(i1)
  // C      dm2=e(i2)
  dm3 = 0.938f;
  dm4 = 0.938f;
  // C     10/24/02 initialize n to 0:
  n = 0;
  // C
  // Cbz3/11/99 neutralk
  // C        if(nchrg.eq.-2.or.nchrg.ge.3) dm3=1.232
  // C        if(nchrg.eq.4) dm4=1.232
  if (nchrg <= -1 || nchrg >= 3) {
    dm3 = 1.232f;
  }
  if (nchrg == -2 || nchrg == 4) {
    dm4 = 1.232f;
  }
  // Cbz3/11/99 neutralk end
  iblock = 0;
  fstate(cmn, iseed, srt, dm3, dm4, px, py, pz, iflag);
  if (iflag < 0) {
    // C           write(60,*)'------------final state fail-------',n
    // C     no anti-kaon production
    ictrl = -1;
    n++;
    return;
  }
  iblock = 12;
  // C Rotate the momenta of particles in the cms of I1 & I2
  // C px(1), py(1), pz(1): momentum of I1
  // C px(2), py(2), pz(2): momentum of I2
  // C px(3), py(3), pz(3): momentum of anti-kaon
  // C px(4), py(4), pz(4): momentum of kaon
  // C
  // C     10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = px(1);
  pyrota = py(1);
  pzrota = pz(1);
  // C        call rotate(pcx,pcy,pcz,px(1),py(1),pz(1))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  px(1) = pxrota;
  py(1) = pyrota;
  pz(1) = pzrota;
  // C
  pxrota = px(2);
  pyrota = py(2);
  pzrota = pz(2);
  // C        call rotate(pcx,pcy,pcz,px(2),py(2),pz(2))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  px(2) = pxrota;
  py(2) = pyrota;
  pz(2) = pzrota;
  // C
  pxrota = px(3);
  pyrota = py(3);
  pzrota = pz(3);
  // C        call rotate(pcx,pcy,pcz,px(3),py(3),pz(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  px(3) = pxrota;
  py(3) = pyrota;
  pz(3) = pzrota;
  // C
  pxrota = px(4);
  pyrota = py(4);
  pzrota = pz(4);
  // C        call rotate(pcx,pcy,pcz,px(4),py(4),pz(4))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  px(4) = pxrota;
  py(4) = pyrota;
  pz(4) = pzrota;
  // C
  nnn += 2;
  // C     K+
  lpion(nnn, irun) = 23;
  if (nchrg == -1 || nchrg == -2) {
    // C        To keep charge conservation. D-n->nnK0K-, D-D- -> nD-K0K-
    // C
    // Cbz3/7/99 neutralk
    // C           lpion(nnn,irun)=24 ! K0
    // Cbz3/7/99 neutralk end
    // C
  }
  // C     aka: rest mass of K
  const float aka = 0.498f;
  epion(nnn, irun) = aka;
  // C     K-
  lpion(nnn - 1, irun) = 21;
  // C     aka: rest mass of K
  epion(nnn - 1, irun) = aka;
  // C Find the momenta of particles in the final state in the nucleus_nucleus
  // C cms frame.   Lorentz transformation into lab frame.
  e1cm = fem::sqrt(fem::pow2(dm3) + fem::pow2(px(1)) + fem::pow2(py(1)) +
                   fem::pow2(pz(1)));
  p1beta = px(1) * betax + py(1) * betay + pz(1) * betaz;
  transf = gamma * (gamma * p1beta / (gamma + 1) + e1cm);
  pt1i1 = betax * transf + px(1);
  pt2i1 = betay * transf + py(1);
  pt3i1 = betaz * transf + pz(1);
  eti1 = dm3;
  // C        lb1   = lb(i1)
  lb1 = 2;
  if (nchrg >= -2 && nchrg <= 1) {
    lb1 = 2;
  }
  // C
  // Cbz3/7/99 neutralk
  if (nchrg == -2 || nchrg == -1) {
    lb1 = 6;
  }
  // Cbz3/7/99 neutralk end
  // C
  // Cbz3/11/99 neutralk
  // C        if(nchrg.eq.2.or.nchrg.eq.3) lb1=1
  // C        if(nchrg.eq.4) lb1=9
  if (nchrg == 1 || nchrg == 2) {
    lb1 = 1;
  }
  if (nchrg == 3 || nchrg == 4) {
    lb1 = 9;
  }
  // Cbz3/11/99 neutralk end
  // C
  // C For second nulceon, same
  e2cm = fem::sqrt(fem::pow2(dm4) + fem::pow2(px(2)) + fem::pow2(py(2)) +
                   fem::pow2(pz(2)));
  p2beta = px(2) * betax + py(2) * betay + pz(2) * betaz;
  transf = gamma * (gamma * p2beta / (gamma + 1) + e2cm);
  pt1i2 = betax * transf + px(2);
  pt2i2 = betay * transf + py(2);
  pt3i2 = betaz * transf + pz(2);
  eti2 = dm4;
  // C        lb2   = lb(i2)
  lb2 = 2;
  // C
  // Cbz3/11/99 neutralk
  // C        if(nchrg.eq.-1.or.nchrg.eq.0) lb2=2
  // C        if(nchrg.eq. 2.or.nchrg.eq.1) lb2=1
  // C        if(nchrg.eq. 4.or.nchrg.eq.3) lb2=9
  // C        if(nchrg.eq.-2) lb2=6
  if (nchrg >= -1 || nchrg <= 1) {
    lb2 = 2;
  }
  if (nchrg == 2 || nchrg == 3) {
    lb2 = 1;
  }
  if (nchrg == 4) {
    lb2 = 9;
  }
  if (nchrg == -2) {
    lb2 = 6;
  }
  // Cbz3/11/99 neutralk end
  // C
  // C        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
  p(1, i1) = pt1i1;
  p(2, i1) = pt2i1;
  p(3, i1) = pt3i1;
  e(i1) = eti1;
  lb(i1) = lb1;
  p(1, i2) = pt1i2;
  p(2, i2) = pt2i2;
  p(3, i2) = pt3i2;
  e(i2) = eti2;
  lb(i2) = lb2;
  // C
  // C                px1 = p(1,i1)
  // C                py1 = p(2,i1)
  // C                pz1 = p(3,i1)
  // C                em1 = e(i1)
  // C                id(i1) = 2
  // C                id(i2) = 2
  // C                id1 = id(i1)
  // C                iblock = 101  ! K(+)K(-) production
  // C Get anti-kaons' momenta and coordinates in nucleus-nucleus cms. frame.
  epcmk = fem::sqrt(fem::pow2(epion(nnn - 1, irun)) + fem::pow2(px(3)) +
                    fem::pow2(py(3)) + fem::pow2(pz(3)));
  betak = px(3) * betax + py(3) * betay + pz(3) * betaz;
  transf = gamma * (gamma * betak / (gamma + 1.f) + epcmk);
  ppion(1, nnn - 1, irun) = betax * transf + px(3);
  ppion(2, nnn - 1, irun) = betay * transf + py(3);
  ppion(3, nnn - 1, irun) = betaz * transf + pz(3);
  rpion(1, nnn - 1, irun) = r(1, i1);
  rpion(2, nnn - 1, irun) = r(2, i1);
  rpion(3, nnn - 1, irun) = r(3, i1);
  // Clin-5/2008:
  dppion(nnn - 1, irun) = dpertp(i1) * dpertp(i2);
  // C Same thing for kaon **************************************
  epcmak = fem::sqrt(fem::pow2(epion(nnn, irun)) + fem::pow2(px(4)) +
                     fem::pow2(py(4)) + fem::pow2(pz(4)));
  betaak = px(4) * betax + py(4) * betay + pz(4) * betaz;
  transf = gamma * (gamma * betaak / (gamma + 1.f) + epcmak);
  ppion(1, nnn, irun) = betax * transf + px(4);
  ppion(2, nnn, irun) = betay * transf + py(4);
  ppion(3, nnn, irun) = betaz * transf + pz(4);
  rpion(1, nnn, irun) = r(1, i2);
  rpion(2, nnn, irun) = r(2, i2);
  rpion(3, nnn, irun) = r(3, i2);
  // Clin-5/2008:
  dppion(nnn, irun) = dpertp(i1) * dpertp(i2);
}

struct lorntz_save {
  float bb;
  float deno3;
  float ga;
  float gam;
  int i;
  float pib;
  float pjb;

  lorntz_save()
      : bb(fem::float0),
        deno3(fem::float0),
        ga(fem::float0),
        gam(fem::float0),
        i(fem::int0),
        pib(fem::float0),
        pjb(fem::float0) {}
};

void lorntz(common& cmn, int const& ilo, arr_cref<float> b, arr_ref<float> pi,
            arr_ref<float> pj) {
  FEM_CMN_SVE(lorntz);
  b(dimension(3));
  pi(dimension(4));
  pj(dimension(4));
  float& bb = sve.bb;
  float& deno3 = sve.deno3;
  float& ga = sve.ga;
  float& gam = sve.gam;
  int& i = sve.i;
  float& pib = sve.pib;
  float& pjb = sve.pjb;
  // C       It uses to perform Lorentz (or inverse Lorentz) transformation
  // C       dimension db(3)
  bb = b(1) * b(1) + b(2) * b(2) + b(3) * b(3);
  deno3 = fem::sqrt(1.f - bb);
  if (deno3 == 0.f) {
    deno3 = 1.e-10f;
  }
  gam = 1.f / deno3;
  ga = gam * gam / (gam + 1.f);
  if (ilo == 1) {
    goto statement_100;
  }
  // C       Lorentz transformation
  pib = pi(1) * b(1) + pi(2) * b(2) + pi(3) * b(3);
  pjb = pj(1) * b(1) + pj(2) * b(2) + pj(3) * b(3);
  // C       drb=drd(1)*b(1)+drd(2)*b(2)+drd(3)*b(3)
  // C       drdb=db(1)*b(1)+db(2)*b(2)+db(3)*b(3)
  FEM_DO_SAFE(i, 1, 3) {
    pi(i) += b(i) * (ga * pib - gam * pi(4));
    pj(i) += b(i) * (ga * pjb - gam * pj(4));
    // C       drd(i)=drd(i)+b(i)*ga*drb
    // C       db(i)=db(i)+b(i)*ga*drdb
  }
  pi(4) = gam * (pi(4) - pib);
  pj(4) = gam * (pj(4) - pjb);
  return;
statement_100:
  // C       inverse Lorentz transformation
  pib = pi(1) * b(1) + pi(2) * b(2) + pi(3) * b(3);
  pjb = pj(1) * b(1) + pj(2) * b(2) + pj(3) * b(3);
  FEM_DO_SAFE(i, 1, 3) {
    pi(i) += b(i) * (ga * pib + gam * pi(4));
    pj(i) += b(i) * (ga * pjb + gam * pj(4));
  }
  pi(4) = gam * (pi(4) + pib);
  pj(4) = gam * (pj(4) + pjb);
}

struct npik_save {
  arr<float> bb;
  float betak;
  float css;
  float e1cm;
  float e2cm;
  float eip;
  float epcmk;
  float eti1;
  float eti2;
  float fai;
  int i;
  int ilo;
  int k;
  int k1;
  int k2;
  int lb1;
  int lb2;
  arr<float> p1;
  float p1beta;
  arr<float> p2;
  float p2beta;
  arr<float> p3;
  float pk;
  float pkmax;
  float pt1i1;
  float pt1i2;
  float pt2i1;
  float pt2i2;
  float pt3i1;
  float pt3i2;
  arr<float> px;
  float px1cm;
  float pxrota;
  arr<float> py;
  float py1cm;
  float pyrota;
  arr<float> pz;
  float pz1cm;
  float pznp;
  float pzrota;
  float rmnp;
  float sss;
  float transf;

  npik_save()
      : bb(dimension(3), fem::fill0),
        betak(fem::float0),
        css(fem::float0),
        e1cm(fem::float0),
        e2cm(fem::float0),
        eip(fem::float0),
        epcmk(fem::float0),
        eti1(fem::float0),
        eti2(fem::float0),
        fai(fem::float0),
        i(fem::int0),
        ilo(fem::int0),
        k(fem::int0),
        k1(fem::int0),
        k2(fem::int0),
        lb1(fem::int0),
        lb2(fem::int0),
        p1(dimension(4), fem::fill0),
        p1beta(fem::float0),
        p2(dimension(4), fem::fill0),
        p2beta(fem::float0),
        p3(dimension(4), fem::fill0),
        pk(fem::float0),
        pkmax(fem::float0),
        pt1i1(fem::float0),
        pt1i2(fem::float0),
        pt2i1(fem::float0),
        pt2i2(fem::float0),
        pt3i1(fem::float0),
        pt3i2(fem::float0),
        px(dimension(4), fem::fill0),
        px1cm(fem::float0),
        pxrota(fem::float0),
        py(dimension(4), fem::fill0),
        py1cm(fem::float0),
        pyrota(fem::float0),
        pz(dimension(4), fem::fill0),
        pz1cm(fem::float0),
        pznp(fem::float0),
        pzrota(fem::float0),
        rmnp(fem::float0),
        sss(fem::float0),
        transf(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....extracted from G. Song's ART expasion including K- interactions
// C.....file `NPIK.FOR'
// C
// C***************************************
// C        subroutine npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
// C     &                  pcx,pcy,pcz,nchrg,ratiok)
void npik(common& cmn, int const& irun, int const& /* iseed */,
          float const& /* dt */, int const& /* nt */, int& ictrl, int const& i1,
          int const& i2, float const& srt, float const& pcx, float const& pcy,
          float const& pcz, int const& /* nchrg */, float const& /* ratiok */,
          int& iblock) {
  FEM_CMN_SVE(npik);
  // COMMON aa
  const int maxstr = 150001;
  arr_cref<float, 2> r(cmn.r, dimension(3, maxstr));
  // COMMON bb
  arr_ref<float, 2> p(cmn.p, dimension(3, maxstr));
  // COMMON cc
  arr_ref<float> e(cmn.e, dimension(maxstr));
  // COMMON ee
  arr_ref<int> lb(cmn.lb, dimension(maxstr));
  // COMMON bg
  float& betax = cmn.betax;
  float& betay = cmn.betay;
  float& betaz = cmn.betaz;
  float& gamma = cmn.gamma;
  // COMMON nn
  int& nnn = cmn.nnn;
  // COMMON pa
  const int maxr = 1;
  arr_ref<float, 3> rpion(cmn.rpion, dimension(3, maxstr, maxr));
  // COMMON pb
  arr_ref<float, 3> ppion(cmn.ppion, dimension(3, maxstr, maxr));
  // COMMON pc
  arr_ref<float, 2> epion(cmn.epion, dimension(maxstr, maxr));
  // COMMON pd
  arr_ref<int, 2> lpion(cmn.lpion, dimension(maxstr, maxr));
  // COMMON rndf77
  int& nseed = cmn.nseed;
  // COMMON dpert
  arr_cref<float> dpertp(cmn.dpertp, dimension(maxstr));
  arr_ref<float, 2> dppion(cmn.dppion, dimension(maxstr, maxr));
  //
  // SAVE
  arr_ref<float> bb(sve.bb, dimension(3));
  float& betak = sve.betak;
  float& css = sve.css;
  float& e1cm = sve.e1cm;
  float& e2cm = sve.e2cm;
  float& eip = sve.eip;
  float& epcmk = sve.epcmk;
  float& eti1 = sve.eti1;
  float& eti2 = sve.eti2;
  float& fai = sve.fai;
  int& i = sve.i;
  int& ilo = sve.ilo;
  int& k = sve.k;
  int& k1 = sve.k1;
  int& k2 = sve.k2;
  int& lb1 = sve.lb1;
  int& lb2 = sve.lb2;
  arr_ref<float> p1(sve.p1, dimension(4));
  float& p1beta = sve.p1beta;
  arr_ref<float> p2(sve.p2, dimension(4));
  float& p2beta = sve.p2beta;
  arr_ref<float> p3(sve.p3, dimension(4));
  float& pk = sve.pk;
  float& pkmax = sve.pkmax;
  float& pt1i1 = sve.pt1i1;
  float& pt1i2 = sve.pt1i2;
  float& pt2i1 = sve.pt2i1;
  float& pt2i2 = sve.pt2i2;
  float& pt3i1 = sve.pt3i1;
  float& pt3i2 = sve.pt3i2;
  float& pxrota = sve.pxrota;
  float& pyrota = sve.pyrota;
  float& pznp = sve.pznp;
  float& pzrota = sve.pzrota;
  float& rmnp = sve.rmnp;
  float& sss = sve.sss;
  float& transf = sve.transf;
  //
  // C
  // C Process: PI + N -> K(-) + ANYTHING
  // C 1.  PI- + P -> P + K0 + K-
  // C 2.  PI+ + N -> P + K+ + K-
  // C 3.  PI0 + P -> P + K+ + K-
  // C 4.  PI0 + N -> P + K0 + K-
  // C 5.  PI0 + N -> N + K+ + K-
  // C 6.  PI- + P -> N + K+ + K-
  // C 7.  PI- + N -> N + K0 + K-
  // C NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
  // C***************************************
  // Cc      SAVE /AA/
  // Cc      SAVE /BB/
  // Cc      SAVE /CC/
  // Cc      SAVE /EE/
  // Cc      SAVE /BG/
  // Cc      SAVE /NN/
  // Cc      SAVE /RUN/
  // Cc      SAVE /PA/
  // Cc      SAVE /PB/
  // Cc      SAVE /PC/
  // Cc      SAVE /PD/
  // Cc      SAVE /RNDF77/
  sve.px1cm = pcx;
  sve.py1cm = pcy;
  sve.pz1cm = pcz;
  ictrl = 1;
  lb1 = lb(i1);
  lb2 = lb(i2);
  k1 = i1;
  k2 = i2;
  // C        k1 must be bayron. k2 be meson. If not, exchange.
  if (lb2 == 1 || lb2 == 2 || (lb2 >= 6 && lb2 <= 13)) {
    k1 = i2;
    k2 = i1;
  }
  // Cbz3/8/99 neutralk
  // Cbz10/12/99
  // C        LB(I1) = 1 + 2 * RANART(NSEED)
  // C        LB(I2) = 23
  lb(k1) = 1 + fem::fint(2 * ranart(nseed));
  lb(k2) = 23;
  // C       pkmax=sqrt((srt**2-(aka+0.938+aka)**2)*(srt**2-(aka+0.938-aka)**2))
  // C     &           /2./srt
  const float aka = 0.498f;
  pkmax = fem::sqrt((fem::pow2(srt) - fem::pow2((aka + 0.938f + aka))) *
                    (fem::pow2(srt) - fem::pow2((aka + 0.938f - aka)))) /
          2.f / srt;
  pk = ranart(nseed) * pkmax;
  // C-----------------------------------------------------
  css = 1.f - 2.f * ranart(nseed);
  sss = fem::sqrt(1.f - fem::pow2(css));
  fai = 2 * 3.1415926f * ranart(nseed);
  p3(1) = pk * sss * fem::cos(fai);
  p3(2) = pk * sss * fem::sin(fai);
  p3(3) = pk * css;
  eip = srt - fem::sqrt(fem::pow2(aka) + fem::pow2(pk));
  rmnp = fem::sqrt(fem::pow2(eip) - fem::pow2(pk));
  FEM_DO_SAFE(i, 1, 3) { bb(i) = -1.f * p3(i) / eip; }
  // C        bb: velocity of the other two particles as a whole.
  pznp = fem::sqrt((fem::pow2(rmnp) - fem::pow2((aka + 0.938f))) *
                   (fem::pow2(rmnp) - fem::pow2((0.938f - aka)))) /
         2.f / rmnp;
  // C-----------------------------------------------------
  css = 1.f - 2.f * ranart(nseed);
  sss = fem::sqrt(1.f - fem::pow2(css));
  fai = 2 * 3.1415926f * ranart(nseed);
  p1(1) = pznp * sss * fem::cos(fai);
  p1(2) = pznp * sss * fem::sin(fai);
  p1(3) = pznp * css;
  p1(4) = fem::sqrt(fem::pow2(0.938f) + fem::pow2(pznp));
  p2(4) = fem::sqrt(fem::pow2(aka) + fem::pow2(pznp));
  FEM_DO_SAFE(i, 1, 3) { p2(i) = -1.f * p1(i); }
  // C        p1,p2: the momenta of the two particles in their cms
  // C        p1: momentum of N or P
  // C        p2: momentum of anti_kaon
  // C        p3: momentum of K0 or K+
  ilo = 1;
  // C        write(61,*)'--------p1,p2',p1,p2
  // C        write(61,*)'--------bb',bb
  lorntz(cmn, ilo, bb, p1, p2);
  // C******* Checking *************
  // C        pxsum = p1(1)+p2(1)+p3(1)
  // C        pysum = p1(2)+p2(2)+p3(2)
  // C        pzsum = p1(3)+p2(3)+p3(3)
  // C        pesum = p1(4)+p2(4)+sqrt(p3(1)**2+p3(2)**2+p3(3)**2+aka**2)-srt
  // C        write(61,*)'---p1,pxsum',p1,pxsum
  // C        write(61,*)'---p2,pysum',p2,pysum
  // C        write(61,*)'---p3,pzsum',p3,pzsum
  // C        write(61,*)'---pesum',pesum
  // C***********************************
  // C
  // C Rotate the momenta of particles in the cms of I1 & I2
  // C px(1), py(1), pz(1): momentum of I1
  // C px(2), py(2), pz(2): momentum of I2
  // C px(3), py(3), pz(3): momentum of anti-kaon
  // C
  // C     10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = p1(1);
  pyrota = p1(2);
  pzrota = p1(3);
  // C        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p1(1) = pxrota;
  p1(2) = pyrota;
  p1(3) = pzrota;
  // C
  pxrota = p2(1);
  pyrota = p2(2);
  pzrota = p2(3);
  // C        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p2(1) = pxrota;
  p2(2) = pyrota;
  p2(3) = pzrota;
  // C
  pxrota = p3(1);
  pyrota = p3(2);
  pzrota = p3(3);
  // C        call rotate(pcx,pcy,pcz,p3(1),p3(2),p3(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p3(1) = pxrota;
  p3(2) = pyrota;
  p3(3) = pzrota;
  // C
  nnn++;
  // C     K(-)
  lpion(nnn, irun) = 21;
  // C     aka: rest mass of K
  epion(nnn, irun) = aka;
  // C Find the momenta of particles in the final state in the nucleus_nucleus
  // C cms frame.   Lorentz transformation into lab frame.
  e1cm = fem::sqrt(fem::pow2(0.938f) + fem::pow2(p1(1)) + fem::pow2(p1(2)) +
                   fem::pow2(p1(3)));
  p1beta = p1(1) * betax + p1(2) * betay + p1(3) * betaz;
  transf = gamma * (gamma * p1beta / (gamma + 1) + e1cm);
  pt1i1 = betax * transf + p1(1);
  pt2i1 = betay * transf + p1(2);
  pt3i1 = betaz * transf + p1(3);
  eti1 = 0.938f;
  lb1 = lb(k1);
  // C
  // C For second nulceon, same
  e2cm = fem::sqrt(fem::pow2(aka) + fem::pow2(p3(1)) + fem::pow2(p3(2)) +
                   fem::pow2(p3(3)));
  p2beta = p3(1) * betax + p3(2) * betay + p3(3) * betaz;
  transf = gamma * (gamma * p2beta / (gamma + 1) + e2cm);
  pt1i2 = betax * transf + p3(1);
  pt2i2 = betay * transf + p3(2);
  pt3i2 = betaz * transf + p3(3);
  eti2 = aka;
  lb2 = lb(k2);
  // C
  // C        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
  // C       k1 stand for nucleon, k2 stand for kaon. lpion stand for Kbar.
  p(1, k1) = pt1i1;
  p(2, k1) = pt2i1;
  p(3, k1) = pt3i1;
  e(k1) = eti1;
  lb(k1) = lb1;
  p(1, k2) = pt1i2;
  p(2, k2) = pt2i2;
  p(3, k2) = pt3i2;
  e(k2) = eti2;
  lb(k2) = lb2;
  // C
  // C                px1 = p(1,i1)
  // C                py1 = p(2,i1)
  // C                pz1 = p(3,i1)
  // C                em1 = e(i1)
  // C                id(i1) = 2
  // C                id(i2) = 2
  // C                id1 = id(i1)
  // C     K(+)K(-) production
  iblock = 101;
  // C Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
  // C  p2:  momentum of anti-kaon.
  // C        epcmk = sqrt(epion(nnn,irun)**2 + p2(1)**2 + p2(2)**2 + p2(3)**2)
  epcmk = fem::sqrt(fem::pow2(epion(nnn, irun)) + fem::pow2(p2(1)) +
                    fem::pow2(p2(2)) + fem::pow2(p2(3)));
  betak = p2(1) * betax + p2(2) * betay + p2(3) * betaz;
  transf = gamma * (gamma * betak / (gamma + 1.f) + epcmk);
  ppion(1, nnn, irun) = betax * transf + p2(1);
  ppion(2, nnn, irun) = betay * transf + p2(2);
  ppion(3, nnn, irun) = betaz * transf + p2(3);
  // Clin-5/2008:
  dppion(nnn, irun) = dpertp(i1) * dpertp(i2);
  // Cbz3/2/99
  // C        write(400,*)'2 ', ppion(1,nnn,irun), ppion(2,nnn,irun),
  // C     &                    ppion(3,nnn,irun), dt*nt, srt
  // Cbz3/2/99end
  // C        write(420,*)ppion(1,nnn,irun), ppion(2,nnn,irun),
  // C     &                    ppion(3,nnn,irun), dt*nt, srt
  k = i2;
  if (lb(i1) == 1 || lb(i1) == 2) {
    k = i1;
  }
  rpion(1, nnn, irun) = r(1, k);
  rpion(2, nnn, irun) = r(2, k);
  rpion(3, nnn, irun) = r(3, k);
}

struct pihypn_save {
  float css;
  float dm3;
  float dm4;
  float e1cm;
  float e2cm;
  float eti1;
  float eti2;
  float fai;
  int i;
  int k1;
  int k2;
  int lb1;
  int lb2;
  arr<float> p1;
  float p1beta;
  arr<float> p2;
  float p2beta;
  float pk;
  float pkmax;
  float pt1i1;
  float pt1i2;
  float pt2i1;
  float pt2i2;
  float pt3i1;
  float pt3i2;
  float px1cm;
  float pxrota;
  float py1cm;
  float pyrota;
  float pz1cm;
  float pzrota;
  float sss;
  float transf;

  pihypn_save()
      : css(fem::float0),
        dm3(fem::float0),
        dm4(fem::float0),
        e1cm(fem::float0),
        e2cm(fem::float0),
        eti1(fem::float0),
        eti2(fem::float0),
        fai(fem::float0),
        i(fem::int0),
        k1(fem::int0),
        k2(fem::int0),
        lb1(fem::int0),
        lb2(fem::int0),
        p1(dimension(4), fem::fill0),
        p1beta(fem::float0),
        p2(dimension(4), fem::fill0),
        p2beta(fem::float0),
        pk(fem::float0),
        pkmax(fem::float0),
        pt1i1(fem::float0),
        pt1i2(fem::float0),
        pt2i1(fem::float0),
        pt2i2(fem::float0),
        pt3i1(fem::float0),
        pt3i2(fem::float0),
        px1cm(fem::float0),
        pxrota(fem::float0),
        py1cm(fem::float0),
        pyrota(fem::float0),
        pz1cm(fem::float0),
        pzrota(fem::float0),
        sss(fem::float0),
        transf(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....extracted from G. Song's ART expasion including K- interactions
// C.....file `PIHYPN.FOR'
// C
// C*****************************************
void pihypn(common& cmn, int const& ielstc, int const& /* irun */,
            int const& /* iseed */, float const& /* dt */, int const& /* nt */,
            int& ictrl, int const& i1, int const& i2, float const& srt,
            float const& pcx, float const& pcy, float const& pcz,
            int const& nchrg, int& iblock) {
  FEM_CMN_SVE(pihypn);
  // COMMON bb
  const int maxstr = 150001;
  arr_ref<float, 2> p(cmn.p, dimension(3, maxstr));
  // COMMON cc
  arr_ref<float> e(cmn.e, dimension(maxstr));
  // COMMON ee
  arr_ref<int> lb(cmn.lb, dimension(maxstr));
  // COMMON bg
  float& betax = cmn.betax;
  float& betay = cmn.betay;
  float& betaz = cmn.betaz;
  float& gamma = cmn.gamma;
  // COMMON rndf77
  int& nseed = cmn.nseed;
  //
  // SAVE
  float& css = sve.css;
  float& dm3 = sve.dm3;
  float& dm4 = sve.dm4;
  float& e1cm = sve.e1cm;
  float& e2cm = sve.e2cm;
  float& eti1 = sve.eti1;
  float& eti2 = sve.eti2;
  float& fai = sve.fai;
  int& i = sve.i;
  int& k1 = sve.k1;
  int& k2 = sve.k2;
  int& lb1 = sve.lb1;
  int& lb2 = sve.lb2;
  arr_ref<float> p1(sve.p1, dimension(4));
  float& p1beta = sve.p1beta;
  arr_ref<float> p2(sve.p2, dimension(4));
  float& p2beta = sve.p2beta;
  float& pk = sve.pk;
  float& pkmax = sve.pkmax;
  float& pt1i1 = sve.pt1i1;
  float& pt1i2 = sve.pt1i2;
  float& pt2i1 = sve.pt2i1;
  float& pt2i2 = sve.pt2i2;
  float& pt3i1 = sve.pt3i1;
  float& pt3i2 = sve.pt3i2;
  float& pxrota = sve.pxrota;
  float& pyrota = sve.pyrota;
  float& pzrota = sve.pzrota;
  float& sss = sve.sss;
  float& transf = sve.transf;
  //
  // C
  // C Process: PI + sigma(or Lambda) -> Kbar + N
  // C NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
  // C*****************************************
  // C
  // C NOTE: for PI + Hyperon: the produced kaons have mass 0.498
  // Cc      SAVE /AA/
  // Cc      SAVE /BB/
  // Cc      SAVE /CC/
  // Cc      SAVE /EE/
  // Cc      SAVE /BG/
  // Cc      SAVE /NN/
  // Cc      SAVE /RUN/
  // Cc      SAVE /PA/
  // Cc      SAVE /PB/
  // Cc      SAVE /PC/
  // Cc      SAVE /PD/
  // Cc      SAVE /RNDF77/
  sve.px1cm = pcx;
  sve.py1cm = pcy;
  sve.pz1cm = pcz;
  ictrl = 1;
  // Csp06/07/01
  const float aka = 0.498f;
  if (ielstc == 1) {
    // C    L/Si + meson -> L/Si + meson
    k1 = i1;
    k2 = i2;
    dm3 = e(k1);
    dm4 = e(k2);
    iblock = 10;
  } else {
    iblock = 12;
    // Csp06/07/01 end
    // C        PI + Sigma(or Lambda) -> Kbar + N
    k1 = i1;
    k2 = i2;
    // C        k1 must be bayron! So if I1 is PI, exchange k1 & k2.
    if (lb(i1) < 14 || lb(i1) > 17) {
      k1 = i2;
      k2 = i1;
    }
    // Cbz3/8/99 neutralk
    lb(k1) = 1 + fem::fint(2 * ranart(nseed));
    if (nchrg == -2) {
      lb(k1) = 6;
    }
    // C     if(nchrg.eq.-1) lb(k1)=2
    // C     if(nchrg.eq. 0) lb(k1)=1
    // C     if(nchrg.eq. 1) lb(k1)=9
    if (nchrg == 2) {
      lb(k1) = 9;
    }
    // Cbz3/8/99 neutralk end
    // C
    // C     K-
    lb(k2) = 21;
    dm3 = 0.938f;
    if (nchrg == -2 || nchrg == 1) {
      dm3 = 1.232f;
    }
    dm4 = aka;
    // C        dm3,dm4: the mass of final state particles.
  }
  // C
  // C*******Now, antikaon will be created.
  // C        call antikaon_fstate(iseed,srt,dm1,dm2,dm3,dm4,px,py,pz,icou1)
  // C        pkmax: the maximum momentum of anti-kaon
  pkmax = fem::sqrt((fem::pow2(srt) - fem::pow2((dm3 + dm4))) *
                    (fem::pow2(srt) - fem::pow2((dm3 - dm4)))) /
          2.f / srt;
  pk = pkmax;
  // C-----------------------------------------------------
  css = 1.f - 2.f * ranart(nseed);
  sss = fem::sqrt(1.f - fem::pow2(css));
  fai = 2 * 3.1415926f * ranart(nseed);
  p1(1) = pk * sss * fem::cos(fai);
  p1(2) = pk * sss * fem::sin(fai);
  p1(3) = pk * css;
  FEM_DO_SAFE(i, 1, 3) { p2(i) = -1.f * p1(i); }
  // C        p1,p2: the momenta of the two particles in their cms
  // C        p1: momentum of kaon
  // C        p2: momentum of Kbar
  // C
  // C Rotate the momenta of particles in the cms of I1 & I2
  // Clin-10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = p1(1);
  pyrota = p1(2);
  pzrota = p1(3);
  // C        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p1(1) = pxrota;
  p1(2) = pyrota;
  p1(3) = pzrota;
  // C
  pxrota = p2(1);
  pyrota = p2(2);
  pzrota = p2(3);
  // C        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p2(1) = pxrota;
  p2(2) = pyrota;
  p2(3) = pzrota;
  // Clin-10/28/02-end
  // C
  // C Find the momenta of particles in the final state in the nucleus_nucleus
  // C cms frame.   Lorentz transformation into lab frame.
  e1cm = fem::sqrt(fem::pow2(dm3) + fem::pow2(p1(1)) + fem::pow2(p1(2)) +
                   fem::pow2(p1(3)));
  p1beta = p1(1) * betax + p1(2) * betay + p1(3) * betaz;
  transf = gamma * (gamma * p1beta / (gamma + 1) + e1cm);
  pt1i1 = betax * transf + p1(1);
  pt2i1 = betay * transf + p1(2);
  pt3i1 = betaz * transf + p1(3);
  eti1 = dm3;
  lb1 = lb(k1);
  // C
  // C For second kaon, same
  e2cm = fem::sqrt(fem::pow2(dm4) + fem::pow2(p2(1)) + fem::pow2(p2(2)) +
                   fem::pow2(p2(3)));
  p2beta = p2(1) * betax + p2(2) * betay + p2(3) * betaz;
  transf = gamma * (gamma * p2beta / (gamma + 1) + e2cm);
  pt1i2 = betax * transf + p2(1);
  pt2i2 = betay * transf + p2(2);
  pt3i2 = betaz * transf + p2(3);
  // Cbz3/2/99
  // C        write(400,*)'3 ', pt1i2, pt2i2, pt3i2, dt*nt, srt
  // Cbz3/2/99end
  // C        write(430,*)pt1i2, pt2i2, pt3i2, dt*nt, srt
  eti2 = dm4;
  lb2 = lb(k2);
  // C
  // C        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
  // C        k1=i1
  // C        k2=i2
  // C       k1 stand for nucleon, k2 stand for kaon.
  p(1, k1) = pt1i1;
  p(2, k1) = pt2i1;
  p(3, k1) = pt3i1;
  e(k1) = eti1;
  lb(k1) = lb1;
  p(1, k2) = pt1i2;
  p(2, k2) = pt2i2;
  p(3, k2) = pt3i2;
  e(k2) = eti2;
  lb(k2) = lb2;
  // C
  // Cc                iblock = 101  ! K(+)K(-) production
  // C Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
}

struct kaonn_save {
  float css;
  float e1cm;
  float e2cm;
  float eee;
  float em1;
  float em2;
  float eti1;
  float eti2;
  float fai;
  int i;
  int k1;
  int k2;
  int lb1;
  int lb2;
  arr<float> p1;
  float p1beta;
  arr<float> p2;
  float p2beta;
  float pk;
  float pkmax;
  float pt1i1;
  float pt1i2;
  float pt2i1;
  float pt2i2;
  float pt3i1;
  float pt3i2;
  float px1cm;
  float pxrota;
  float py1cm;
  float pyrota;
  float pz1cm;
  float pzrota;
  float rrr;
  float sss;
  float transf;

  kaonn_save()
      : css(fem::float0),
        e1cm(fem::float0),
        e2cm(fem::float0),
        eee(fem::float0),
        em1(fem::float0),
        em2(fem::float0),
        eti1(fem::float0),
        eti2(fem::float0),
        fai(fem::float0),
        i(fem::int0),
        k1(fem::int0),
        k2(fem::int0),
        lb1(fem::int0),
        lb2(fem::int0),
        p1(dimension(4), fem::fill0),
        p1beta(fem::float0),
        p2(dimension(4), fem::fill0),
        p2beta(fem::float0),
        pk(fem::float0),
        pkmax(fem::float0),
        pt1i1(fem::float0),
        pt1i2(fem::float0),
        pt2i1(fem::float0),
        pt2i2(fem::float0),
        pt3i1(fem::float0),
        pt3i2(fem::float0),
        px1cm(fem::float0),
        pxrota(fem::float0),
        py1cm(fem::float0),
        pyrota(fem::float0),
        pz1cm(fem::float0),
        pzrota(fem::float0),
        rrr(fem::float0),
        sss(fem::float0),
        transf(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....extracted from G. Song's ART expasion including K- interactions
// C.....file `KAONN.FOR'
// C
// C***************************************
void kaonn(common& cmn, float const& brel, float const& brsgm,
           int const& /* irun */, int const& /* iseed */, float const& /* dt */,
           int const& /* nt */, int& ictrl, int const& i1, int const& i2,
           int& iblock, float const& srt, float const& pcx, float const& pcy,
           float const& pcz, int const& /* nchrg */) {
  FEM_CMN_SVE(kaonn);
  // COMMON bb
  const int maxstr = 150001;
  arr_ref<float, 2> p(cmn.p, dimension(3, maxstr));
  // COMMON cc
  arr_ref<float> e(cmn.e, dimension(maxstr));
  // COMMON ee
  arr_ref<int> lb(cmn.lb, dimension(maxstr));
  // COMMON bg
  float& betax = cmn.betax;
  float& betay = cmn.betay;
  float& betaz = cmn.betaz;
  float& gamma = cmn.gamma;
  // COMMON rndf77
  int& nseed = cmn.nseed;
  //
  // SAVE
  float& css = sve.css;
  float& e1cm = sve.e1cm;
  float& e2cm = sve.e2cm;
  float& em1 = sve.em1;
  float& em2 = sve.em2;
  float& eti1 = sve.eti1;
  float& eti2 = sve.eti2;
  float& fai = sve.fai;
  int& i = sve.i;
  int& k1 = sve.k1;
  int& k2 = sve.k2;
  int& lb1 = sve.lb1;
  int& lb2 = sve.lb2;
  arr_ref<float> p1(sve.p1, dimension(4));
  float& p1beta = sve.p1beta;
  arr_ref<float> p2(sve.p2, dimension(4));
  float& p2beta = sve.p2beta;
  float& pk = sve.pk;
  float& pkmax = sve.pkmax;
  float& pt1i1 = sve.pt1i1;
  float& pt1i2 = sve.pt1i2;
  float& pt2i1 = sve.pt2i1;
  float& pt2i2 = sve.pt2i2;
  float& pt3i1 = sve.pt3i1;
  float& pt3i2 = sve.pt3i2;
  float& pxrota = sve.pxrota;
  float& pyrota = sve.pyrota;
  float& pzrota = sve.pzrota;
  float& rrr = sve.rrr;
  float& sss = sve.sss;
  float& transf = sve.transf;
  //
  // C
  // C Process: PI + sigma(or Lambda) <- Kbar + N
  // C NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
  // C***************************************
  // Cc      SAVE /AA/
  // Cc      SAVE /BB/
  // Cc      SAVE /CC/
  // Cc      SAVE /EE/
  // Cc      SAVE /BG/
  // Cc      SAVE /NN/
  // Cc      SAVE /RUN/
  // Cc      SAVE /PA/
  // Cc      SAVE /PB/
  // Cc      SAVE /PC/
  // Cc      SAVE /PD/
  // Cc      SAVE /RNDF77/
  sve.px1cm = pcx;
  sve.py1cm = pcy;
  sve.pz1cm = pcz;
  ictrl = 1;
  // C        ratio: used for isospin decision.
  k1 = i1;
  k2 = i2;
  // C        k1 must be bayron! So if I1 is Kaon, exchange k1 & k2.
  if (e(i1) < 0.5f && e(i1) > 0.01f) {
    k1 = i2;
    k2 = i1;
  }
  // C** note: for print out only *******************************
  // C     record kaon's mass
  sve.eee = e(k2);
  // C** end **************
  rrr = ranart(nseed);
  const float asa = 1.1974f;
  const float ala = 1.1157f;
  if (rrr < brel) {
    // C       Kbar + N -> Kbar + N
    lb1 = lb(k1);
    lb2 = lb(k2);
    em1 = e(k1);
    em2 = e(k2);
    iblock = 10;
  } else {
    iblock = 12;
    if (rrr < (brel + brsgm)) {
      // C        nchrg: Net charges of the two incoming particles.
      // C           Kbar + N -> Sigma + PI
      em1 = asa;
      em2 = 0.138f;
      // C
      // Cbz3/8/99 neutralk
      lb1 = 15 + fem::fint(3 * ranart(nseed));
      lb2 = 3 + fem::fint(3 * ranart(nseed));
    } else {
      // C           Kbar + N -> Lambda + PI
      em1 = ala;
      em2 = 0.138f;
      // C     LAmbda
      lb1 = 14;
      // Cbz3/8/99 neutralk
      lb2 = 3 + fem::fint(3 * ranart(nseed));
      // C           if(nchrg.eq.1)  lb2=5  ! K- + D++ -> Lambda + PI+
      // C           if(nchrg.eq.0)  lb2=4  ! K- + p(D+,N*+) -> Lambda + PI0
      // C          if(nchrg.eq.-1) lb2=3 ! K- + n(D,N*) -> Lambda + PI-
      // Cbz3/8/99 neutralk
      // C
    }
  }
  lb(k1) = lb1;
  lb(k2) = lb2;
  // C
  // C*******Now, antikaon will be created.
  // C        call antikaon_fstate(iseed,srt,dm1,dm2,dm3,dm4,px,py,pz,icou1)
  // C        pkmax: the maximum momentum of anti-kaon
  // C        write(63,*)'srt,em1,em2',srt,em1,em2
  // C        write(63,*)'-srt,em1,em2',srt,em1,em2
  pkmax = fem::sqrt((fem::pow2(srt) - fem::pow2((em1 + em2))) *
                    (fem::pow2(srt) - fem::pow2((em1 - em2)))) /
          2.f / srt;
  pk = pkmax;
  // C-----------------------------------------------------
  css = 1.f - 2.f * ranart(nseed);
  sss = fem::sqrt(1.f - fem::pow2(css));
  fai = 2 * 3.1415926f * ranart(nseed);
  p1(1) = pk * sss * fem::cos(fai);
  p1(2) = pk * sss * fem::sin(fai);
  p1(3) = pk * css;
  FEM_DO_SAFE(i, 1, 3) { p2(i) = -1.f * p1(i); }
  // C        p1,p2: the momenta of the two particles in their cms
  // C        p1: momentum of kaon
  // C        p2: momentum of Kbar
  // C
  // C Rotate the momenta of particles in the cms of I1 & I2
  // C
  // Clin-10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = p1(1);
  pyrota = p1(2);
  pzrota = p1(3);
  // C        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p1(1) = pxrota;
  p1(2) = pyrota;
  p1(3) = pzrota;
  // C
  pxrota = p2(1);
  pyrota = p2(2);
  pzrota = p2(3);
  // C        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
  rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota);
  p2(1) = pxrota;
  p2(2) = pyrota;
  p2(3) = pzrota;
  // Clin-10/28/02-end
  // C
  // C Find the momenta of particles in the final state in the nucleus_nucleus
  // C cms frame.   Lorentz transformation into lab frame.
  e1cm = fem::sqrt(fem::pow2(em1) + fem::pow2(p1(1)) + fem::pow2(p1(2)) +
                   fem::pow2(p1(3)));
  p1beta = p1(1) * betax + p1(2) * betay + p1(3) * betaz;
  transf = gamma * (gamma * p1beta / (gamma + 1) + e1cm);
  pt1i1 = betax * transf + p1(1);
  pt2i1 = betay * transf + p1(2);
  pt3i1 = betaz * transf + p1(3);
  eti1 = em1;
  // C
  // C For second kaon, same
  e2cm = fem::sqrt(fem::pow2(em2) + fem::pow2(p2(1)) + fem::pow2(p2(2)) +
                   fem::pow2(p2(3)));
  p2beta = p2(1) * betax + p2(2) * betay + p2(3) * betaz;
  transf = gamma * (gamma * p2beta / (gamma + 1) + e2cm);
  pt1i2 = betax * transf + p2(1);
  pt2i2 = betay * transf + p2(2);
  pt3i2 = betaz * transf + p2(3);
  eti2 = em2;
  // C
  // C        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
  // C        k1=i1
  // C        k2=i2
  // C       k1 stand for bayron, k2 stand for meson.
  p(1, k1) = pt1i1;
  p(2, k1) = pt2i1;
  p(3, k1) = pt3i1;
  e(k1) = eti1;
  p(1, k2) = pt1i2;
  p(2, k2) = pt2i2;
  p(3, k2) = pt3i2;
  e(k2) = eti2;
  // C
  // Cc                iblock = 101  ! K(+)K(-) production
  // C Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
}

struct newka_save {
  float bmass;
  float brel;
  float brsgm;
  float brsig;
  float ds;
  float dsr;
  float ec;
  float em1;
  float em2;
  int ic;
  int ielstc;
  int ik;
  int ik0;
  int ik1;
  int ik2;
  int ik3;
  int il;
  int im;
  int im3;
  int im4;
  int in;
  int inpion;
  int ipipi;
  int lb1;
  bool lb1bn;
  bool lb1bn0;
  bool lb1bn1;
  bool lb1mn;
  bool lb1mn0;
  bool lb1mn1;
  bool lb1mn2;
  int lb2;
  bool lb2bn;
  bool lb2bn0;
  bool lb2bn1;
  bool lb2mn;
  bool lb2mn0;
  bool lb2mn1;
  bool lb2mn2;
  int nchrg;
  float pkaon;
  float px1cm;
  float py1cm;
  float pz1cm;
  float ratiok;
  float sgsum;
  float sgsum1;
  float sgsum3;
  float sig;
  float sigela;
  float sigma0;
  float sigsgm;

  newka_save()
      : bmass(fem::float0),
        brel(fem::float0),
        brsgm(fem::float0),
        brsig(fem::float0),
        ds(fem::float0),
        dsr(fem::float0),
        ec(fem::float0),
        em1(fem::float0),
        em2(fem::float0),
        ic(fem::int0),
        ielstc(fem::int0),
        ik(fem::int0),
        ik0(fem::int0),
        ik1(fem::int0),
        ik2(fem::int0),
        ik3(fem::int0),
        il(fem::int0),
        im(fem::int0),
        im3(fem::int0),
        im4(fem::int0),
        in(fem::int0),
        inpion(fem::int0),
        ipipi(fem::int0),
        lb1(fem::int0),
        lb1bn(fem::bool0),
        lb1bn0(fem::bool0),
        lb1bn1(fem::bool0),
        lb1mn(fem::bool0),
        lb1mn0(fem::bool0),
        lb1mn1(fem::bool0),
        lb1mn2(fem::bool0),
        lb2(fem::int0),
        lb2bn(fem::bool0),
        lb2bn0(fem::bool0),
        lb2bn1(fem::bool0),
        lb2mn(fem::bool0),
        lb2mn0(fem::bool0),
        lb2mn1(fem::bool0),
        lb2mn2(fem::bool0),
        nchrg(fem::int0),
        pkaon(fem::float0),
        px1cm(fem::float0),
        py1cm(fem::float0),
        pz1cm(fem::float0),
        ratiok(fem::float0),
        sgsum(fem::float0),
        sgsum1(fem::float0),
        sgsum3(fem::float0),
        sig(fem::float0),
        sigela(fem::float0),
        sigma0(fem::float0),
        sigsgm(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....extracted from G. Song's ART expasion including K- interactions
// C.....file `NEWKAON.FOR'
// C
// C     5/01/03 send iblock value into art1f.f, necessary for resonance
// studies: C        subroutine newka(icase,irun,iseed,dt,nt,ictrl,i1,i2, C     &
// srt,pcx,pcy,pcz)
void newka(common& cmn, int& icase, int const& irun, int const& iseed,
           float const& dt, int const& nt, int& ictrl, int const& i1,
           int const& i2, float const& srt, float const& pcx, float const& pcy,
           float const& pcz, int& iblock) {
  FEM_CMN_SVE(newka);
  // COMMON cc
  const int maxstr = 150001;
  arr_cref<float> e(cmn.e, dimension(maxstr));
  // COMMON ee
  arr_cref<int> lb(cmn.lb, dimension(maxstr));
  // COMMON rndf77
  int& nseed = cmn.nseed;
  //
  // SAVE
  float& bmass = sve.bmass;
  float& brel = sve.brel;
  float& brsgm = sve.brsgm;
  float& brsig = sve.brsig;
  float& ds = sve.ds;
  float& dsr = sve.dsr;
  float& ec = sve.ec;
  float& em1 = sve.em1;
  float& em2 = sve.em2;
  int& ic = sve.ic;
  int& ielstc = sve.ielstc;
  int& ik = sve.ik;
  int& ik0 = sve.ik0;
  int& ik1 = sve.ik1;
  int& ik2 = sve.ik2;
  int& ik3 = sve.ik3;
  int& il = sve.il;
  int& im = sve.im;
  int& im3 = sve.im3;
  int& im4 = sve.im4;
  int& in = sve.in;
  int& inpion = sve.inpion;
  int& ipipi = sve.ipipi;
  int& lb1 = sve.lb1;
  bool& lb1bn = sve.lb1bn;
  bool& lb1bn0 = sve.lb1bn0;
  bool& lb1bn1 = sve.lb1bn1;
  bool& lb1mn = sve.lb1mn;
  bool& lb1mn0 = sve.lb1mn0;
  bool& lb1mn1 = sve.lb1mn1;
  bool& lb1mn2 = sve.lb1mn2;
  int& lb2 = sve.lb2;
  bool& lb2bn = sve.lb2bn;
  bool& lb2bn0 = sve.lb2bn0;
  bool& lb2bn1 = sve.lb2bn1;
  bool& lb2mn = sve.lb2mn;
  bool& lb2mn0 = sve.lb2mn0;
  bool& lb2mn1 = sve.lb2mn1;
  bool& lb2mn2 = sve.lb2mn2;
  int& nchrg = sve.nchrg;
  float& pkaon = sve.pkaon;
  float& px1cm = sve.px1cm;
  float& py1cm = sve.py1cm;
  float& pz1cm = sve.pz1cm;
  float& ratiok = sve.ratiok;
  float& sgsum = sve.sgsum;
  float& sgsum1 = sve.sgsum1;
  float& sgsum3 = sve.sgsum3;
  float& sig = sve.sig;
  float& sigela = sve.sigela;
  float& sigma0 = sve.sigma0;
  float& sigsgm = sve.sigsgm;
  //
  // Cc      SAVE /AA/
  // Cc      SAVE /BB/
  // Cc      SAVE /CC/
  // Cc      SAVE /EE/
  // Cc      SAVE /BG/
  // Cc      SAVE /NN/
  // Cc      SAVE /RUN/
  // Cc      SAVE /PA/
  // Cc      SAVE /PB/
  // Cc      SAVE /PC/
  // Cc      SAVE /PD/
  // Cc      SAVE /RNDF77/
  // C
  // Cbz3/7/99 neutralk
  // C        logical lb1bn1, lb2bayon1, lb1bn0, lb2bn0
  // Cbz3/7/99 neutralk end
  icase = -1;
  // C        icase: flag for the type of reaction that is going to happen.
  // C        icase=-1,  no desired reaction, return to main program.
  // C              1,  NN,ND,DD
  // C              2,  PI+N, PI+D
  // C              3,  K(-) absorption.
  nchrg = -100;
  // C        nchrg: Net charges of the two incoming particles.
  ictrl = 1;
  lb1 = lb(i1);
  lb2 = lb(i2);
  em1 = e(i1);
  em2 = e(i2);
  lb1bn = lb1 == 1 || lb1 == 2 || (lb1 > 5 && lb1 <= 13);
  lb2bn = lb2 == 1 || lb2 == 2 || (lb2 > 5 && lb2 <= 13);
  lb1bn0 = lb1 == 2 || lb1 == 7 || lb1 == 10 || lb1 == 12;
  lb2bn0 = lb2 == 2 || lb2 == 7 || lb2 == 10 || lb2 == 12;
  lb1bn1 = lb1 == 1 || lb1 == 8 || lb1 == 11 || lb1 == 13;
  lb2bn1 = lb2 == 1 || lb2 == 8 || lb2 == 11 || lb2 == 13;
  lb1mn = em1 < 0.2f || lb1 == 0 || (lb1 >= 25 && lb1 <= 29);
  lb2mn = em2 < 0.2f || lb2 == 0 || (lb2 >= 25 && lb2 <= 29);
  lb1mn0 = lb1 == 0 || lb1 == 4 || lb1 == 26 || lb1 == 28 || lb1 == 29;
  lb2mn0 = lb2 == 0 || lb2 == 4 || lb2 == 26 || lb2 == 28 || lb2 == 29;
  lb1mn1 = lb1 == 5 || lb1 == 27;
  lb2mn1 = lb2 == 5 || lb2 == 27;
  lb1mn2 = lb1 == 3 || lb1 == 25;
  lb2mn2 = lb2 == 3 || lb2 == 25;
  // C
  // C        1. consider N+N, N+Resonance, R + R reactions
  if (lb1bn && lb2bn) {
    // C     NN,ND,DD:
    icase = 1;
    // C     total cross section
    sig = 40.f;
    if (lb1 == 9 && lb2 == 9) {
      nchrg = 4;
    }
    if ((lb1bn1 && lb2 == 9) || (lb2bn1 && lb1 == 9)) {
      nchrg = 3;
    }
    if ((lb1bn0 && lb2 == 9) || (lb2bn0 && lb1 == 9) || (lb1bn1 && lb2bn1)) {
      nchrg = 2;
    }
    if ((lb1bn1 && lb2bn0) || (lb1 == 6 && lb2 == 9) || (lb2bn1 && lb1bn0) ||
        (lb2 == 6 && lb1 == 9)) {
      nchrg = 1;
    }
    if ((lb1bn0 && lb2bn0) || (lb1bn1 && lb2 == 6) || (lb2bn1 && lb1 == 6)) {
      nchrg = 0;
    }
    if ((lb1bn0 && lb2 == 6) || (lb2bn0 && lb1 == 6)) {
      nchrg = -1;
    }
    if (lb1 == 6 && lb2 == 6) {
      nchrg = -2;
    }
    // C     brsig = x2kaon_no_isospin(srt)
    if (nchrg >= -1 && nchrg <= 2) {
      // C     K,Kbar prduction x sect.
      brsig = x2kaon(cmn, srt);
    } else {
      brsig = 0.0f;
      // C                if(nchrg.eq.-2.or.nchrg.eq.3) then
      // C                   brsig = x2kaon(srt+0.938-1.232)
      // C                else
      // C     nchrg=4
      // C                   brsig = x2kaon(srt+2.*(0.938-1.232))
      // C                endif
    }
    // C
    // Cbz3/7/99 neutralk
    brsig = 2.0f * brsig;
    // Cbz3/7/99 neutralk end
    // C
  }
  // C
  // C        2. consider PI(meson:eta,omega,rho,phi) + N(N*,D)
  if ((lb1bn && lb2mn) || (lb2bn && lb1mn)) {
    // C     PN,PD
    icase = 2;
    sig = 20.f;
    sigma0 = pinsg0(cmn, srt);
    brsig = 0.0f;
    if ((lb1bn1 && lb2mn0) || (lb2bn1 && lb1mn0) || (lb1bn0 && lb2mn1) ||
        (lb2bn0 && lb1mn1) || (lb1 == 9 && lb2mn2) || (lb2 == 9 && lb1mn2)) {
      nchrg = 1;
      // Cbz3/2/99/song
      // C                if(lb1bn1.or.lb2bn1) brsig=2.0*sigma0
      // C                if(lb1bn0.or.lb2bn0) brsig=0.5*sigma0
      if (lb1bn1 || lb2bn1) {
        brsig = 0.5f * sigma0;
      }
      if (lb1bn0 || lb2bn0) {
        brsig = 2.0f * sigma0;
      }
      // Cbz3/2/99/song end
      // C                if(lb1.eq.9.or.lb2.eq.9) brsig=1.5*sigma0
    }
    if ((lb1bn0 && lb2mn0) || (lb2bn0 && lb1mn0) || (lb1bn1 && lb2mn2) ||
        (lb2bn1 && lb1mn2) || (lb1 == 6 && lb2mn1) || (lb2 == 6 && lb1mn1)) {
      nchrg = 0;
      if (lb1bn1 || lb2bn1) {
        // Cbz3/2/99/song
        // C                  brsig=1.5*sigma0
        brsig = 3.0f * sigma0;
        // Cbz3/2/99/song end
        // Cbz3/11/99/song
        // C                  ratiok = 1./3.
        ratiok = 2.f / 3.f;
        // Cbz3/11/99/song end
        // C
        // C                  ratiok: the ratio of channels: ->nK+k- vs. ->
        // pK0K-
      }
      if (lb1bn0 || lb2bn0) {
        brsig = 2.5f * sigma0;
        // Cbz3/2/99/song
        // C                  ratiok = 0.8
        ratiok = 0.2f;
        // Cbz3/2/99/song end
      }
      // C                if(lb1.eq.6.or.lb2.eq.6) then
      // C     lb=6 : D-
      // C                  brsig=1.5*sigma0
      // C                  ratiok = 0.5
      // C                endif
    }
    if ((lb1bn0 && lb2mn2) || (lb2bn0 && lb1mn2) || (lb1 == 6 && lb2mn0) ||
        (lb2 == 6 && lb1mn0)) {
      nchrg = -1;
      if (lb1bn0 || lb2bn0) {
        brsig = sigma0;
      }
      // C                if(lb1.eq.6.or.lb2.eq.6) brsig=sigma0
    }
    // C          if((lb1.eq.6.and.lb2mn2).or.(lb2.eq.6.and.lb1mn2))then
    // C                nchrg=-2
    // C          endif
    // C          if((lb1bn1.and.lb2mn1).or.(lb2bn1.and.lb1mn1)
    // C    &           .or.(lb1.eq.9.and.lb2mn0).or.(lb2.eq.9.and.lb1mn0)) then
    // C                nchrg=2
    // C          endif
    // C
    // Cbz3/11/99 neutralk
    if ((lb1 == 6 && lb2mn2) || (lb2 == 6 && lb1mn2)) {
      nchrg = -2;
    }
    // Cbz3/11/99 neutralk
    // Cbz3/8/99 neutralk
    if ((lb1bn1 && lb2mn1) || (lb2bn1 && lb1mn1) || (lb1 == 9 && lb2mn0) ||
        (lb2 == 9 && lb1mn0)) {
      nchrg = 2;
    }
    // Cbz3/8/99 neutralk end
    // C
    // Cbz3/7/99 neutralk
    if (nchrg >= -2 && nchrg <= 2) {
      brsig = 3.0f * sigma0;
    }
    // Cbz3/7/99 neutralk end
    // C
  }
  // C
  // C        3. consider K- + N(N*,D) absorption.
  // C        if((lb1bn.and.lb2.eq.21).OR.(lb2bn.and.lb1.eq.21)) then
  const float aka = 0.498f;
  if ((lb1bn && (lb2 == 21 || lb2 == -30)) ||
      (lb2bn && (lb1 == 21 || lb1 == -30))) {
    // C          bmass=em1+em2-aka
    bmass = 0.938f;
    if (srt <= (bmass + aka)) {
      // Cbz3/2/99
      // C write(100,*)'--lb1,lb2,em1,em2,srt',lb1,lb2,em1,em2,srt Cbz3/2/99end
      pkaon = 0.f;
    } else {
      pkaon = fem::sqrt(
          fem::pow2(((fem::pow2(srt) - (fem::pow2(aka) + fem::pow2(bmass))) /
                     2.f / bmass)) -
          fem::pow2(aka));
    }
    sig = 0.f;
    if (lb1 == 1 || lb2 == 1 || lb1 == 8 || lb2 == 8 || lb1 == 11 ||
        lb2 == 11 || lb1 == 13 || lb2 == 13) {
      // C          K- + (D+,N*+)p ->
      nchrg = 0;
      sigela = akpel(cmn, pkaon);
      sigsgm = 3.f * akpsgm(cmn, pkaon);
      sig = sigela + sigsgm + akplam(cmn, pkaon);
    }
    if (lb1 == 2 || lb2 == 2 || lb1 == 7 || lb2 == 7 || lb1 == 10 ||
        lb2 == 10 || lb1 == 12 || lb2 == 12) {
      // C          K- + (D0, N*0)n ->
      nchrg = -1;
      sigela = aknel(cmn, pkaon);
      sigsgm = 2.f * aknsgm(cmn, pkaon);
      sig = sigela + sigsgm + aknlam(cmn, pkaon);
    }
    if (lb1 == 6 || lb2 == 6) {
      // C     K- + D-
      nchrg = -2;
      sigela = aknel(cmn, pkaon);
      sigsgm = aknsgm(cmn, pkaon);
      sig = sigela + sigsgm;
    }
    if (lb1 == 9 || lb2 == 9) {
      // C     K- + D++
      nchrg = 1;
      sigela = akpel(cmn, pkaon);
      sigsgm = 2.f * akpsgm(cmn, pkaon);
      sig = sigela + sigsgm + akplam(cmn, pkaon);
    }
    // C
    // Cbz3/8/99 neutralk
    sigela = 0.5f * (akpel(cmn, pkaon) + aknel(cmn, pkaon));
    sigsgm = 1.5f * akpsgm(cmn, pkaon) + aknsgm(cmn, pkaon);
    sig = sigela + sigsgm + akplam(cmn, pkaon);
    // Cbz3/8/99 neutralk end
    // C
    if (sig > 1.e-7f) {
      // C     K(-) + N reactions
      icase = 3;
      brel = sigela / sig;
      brsgm = sigsgm / sig;
      // C              branch_lambda=akNlam(pkaon)/sig
      brsig = sig;
    }
  }
  // C
  // C        4. meson + hyperon -> K- + N
  // C        if(((lb1.ge.14.and.lb1.le.17).and.lb2mn).OR.
  // C     &     ((lb2.ge.14.and.lb2.le.17).and.lb1mn)) then
  if (((lb1 >= 14 && lb1 <= 17) && (lb2 >= 3 && lb2 <= 5)) ||
      ((lb2 >= 14 && lb2 <= 17) && (lb1 >= 3 && lb1 <= 5))) {
    // C        first classify the reactions due to total charge.
    nchrg = -100;
    if ((lb1 == 15 && (lb2 == 3 || lb2 == 25)) ||
        (lb2 == 15 && (lb1 == 3 || lb1 == 25))) {
      nchrg = -2;
      // C     D-
      bmass = 1.232f;
    }
    if ((lb1 == 15 && lb2mn0) || (lb2 == 15 && lb1mn0) ||
        ((lb1 == 14 || lb1 == 16) && (lb2 == 3 || lb2 == 25)) ||
        ((lb2 == 14 || lb2 == 16) && (lb1 == 3 || lb1 == 25))) {
      nchrg = -1;
      // C     n
      bmass = 0.938f;
    }
    if ((lb1 == 15 && (lb2 == 5 || lb2 == 27)) ||
        (lb2 == 15 && (lb1 == 5 || lb1 == 27)) ||
        (lb1 == 17 && (lb2 == 3 || lb2 == 25)) ||
        (lb2 == 17 && (lb1 == 3 || lb1 == 25)) ||
        ((lb1 == 14 || lb1 == 16) && lb2mn0) ||
        ((lb2 == 14 || lb2 == 16) && lb1mn0)) {
      nchrg = 0;
      // C     p
      bmass = 0.938f;
    }
    if ((lb1 == 17 && lb2mn0) || (lb2 == 17 && lb1mn0) ||
        ((lb1 == 14 || lb1 == 16) && (lb2 == 5 || lb2 == 27)) ||
        ((lb2 == 14 || lb2 == 16) && (lb1 == 5 || lb1 == 27))) {
      nchrg = 1;
      // C     D++
      bmass = 1.232f;
    }
    sig = 0.f;
    if (nchrg != -100 && srt > (aka + bmass)) {
      // C     PI+sigma or PI + Lambda => Kbar + N reactions
      icase = 4;
      // C pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
      pkaon = fem::sqrt(
          fem::pow2(((fem::pow2(srt) - (fem::pow2(aka) + fem::pow2(0.938f))) /
                     2.f / 0.938f)) -
          fem::pow2(aka));
      // C     lambda + Pi
      if (lb1 == 14 || lb2 == 14) {
        if (nchrg >= 0) {
          sigma0 = akplam(cmn, pkaon);
        }
        if (nchrg < 0) {
          sigma0 = aknlam(cmn, pkaon);
        }
        // C     sigma + pi
      } else {
        // C     K-p or K-D++
        if (nchrg >= 0) {
          sigma0 = akpsgm(cmn, pkaon);
        }
        // C     K-n or K-D-
        if (nchrg < 0) {
          sigma0 = aknsgm(cmn, pkaon);
        }
        // C
        // Cbz3/8/99 neutralk
        sigma0 = 1.5f * akpsgm(cmn, pkaon) + aknsgm(cmn, pkaon);
        // Cbz3/8/99 neutralk end
        // C
      }
      sig = (fem::pow2(srt) - fem::pow2((aka + bmass))) *
            (fem::pow2(srt) - fem::pow2((aka - bmass))) /
            (fem::pow2(srt) - fem::pow2((em1 + em2))) /
            (fem::pow2(srt) - fem::pow2((em1 - em2))) * sigma0;
      // Cbz3/8/99 neutralk
      // C     if(nchrg.eq.-2.or.nchrg.eq.1) sig=2.*sig K-D++, K-D-
      // C     K0barD++, K-D-
      if (nchrg == -2 || nchrg == 2) {
        sig = 2.f * sig;
      }
      // C
      // Cbz3/8/99 neutralk end
      // C
      // C             the factor 2 comes from spin of delta, which is 3/2
      // C             detailed balance. copy from Page 423 of N.P. A614 1997
      // C
      // Cbz3/8/99 neutralk
      if (lb1 == 14 || lb2 == 14) {
        sig = 4.0f / 3.0f * sig;
      } else if (nchrg == -2 || nchrg == 2) {
        sig = 8.0f / 9.0f * sig;
      } else {
        sig = 4.0f / 9.0f * sig;
      }
      // Cbz3/8/99 neutralk end
      brsig = sig;
      if (sig < 1.e-7f) {
        sig = 1.e-7f;
      }
    }
    // Csp05/07/01
    // C comment icase=4 statement below if only inelastic
    // C     PI+L/Si => Kbar + N  OR ELASTIC SCATTERING
    icase = 4;
    brsig = sig;
    // C     elastic xsecn of 10mb
    sigela = 10.f;
    sig += sigela;
    brel = sigela / sig;
    // Cc          brsig = sig
    // Csp05/07/01 end
  }
  // C
  // C        if(em2.lt.0.2.and.em1.lt.0.2) then
  // C     PI + PI
  // C             icase=5
  // C     assumed PI PI total x section.
  // C              sig=50.
  // C     Mk + Mkbar
  // C              s0=aka+aka
  // C              brsig = 0.
  // C              if(srt.gt.s0) brsig = 2.7*(1.-s0**2/srt**2)**0.76
  // C              x section for PIPI->KKbar   PRC43 (1991) 1881
  // C        endif
  if (icase == -1) {
    ictrl = -1;
    return;
  }
  px1cm = pcx;
  py1cm = pcy;
  pz1cm = pcz;
  ds = fem::sqrt(sig / 31.4f);
  dsr = ds + 0.1f;
  ec = fem::pow2((em1 + em2 + 0.02f));
  // C        ec=3.59709
  // C        if((e(i1).ge.1.).and.(e(i2).ge.1.)) ec = 4.75
  // C
  distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm);
  if (ic == -1) {
    // C     no anti-kaon production
    ictrl = -1;
    // C           in=in+1
    // C           write(60,*)'--------------distance-----',in
    return;
  }
  // C
  // Clin-10/24/02 set to 0: ik,ik0-3,il,im,im3-4,in,inpion,ipipi,
  // C     sgsum,sgsum1,sgsum3:
  ik = 0;
  ik0 = 0;
  ik1 = 0;
  ik2 = 0;
  ik3 = 0;
  il = 0;
  im = 0;
  im3 = 0;
  im4 = 0;
  in = 0;
  inpion = 0;
  ipipi = 0;
  sgsum = 0.f;
  sgsum1 = 0.f;
  sgsum3 = 0.f;
  if (icase == 1) {
    ik++;
    if (srt > 2.8639f) {
      ik0++;
      if (em1 < 1.0f && em2 < 1.0f) {
        ik1++;
        sgsum1 += brsig;
        // C                        ratio_1=sgsum1/ik1/40.
      }
      if (em1 > 1.0f && em2 > 1.0f) {
        ik3++;
        sgsum3 += brsig;
        // C                        ratio_3=sgsum3/ik3/40.
      }
      if (em1 > 1.0f && em2 < 1.0f) {
        ik2++;
      }
      if (em1 < 1.0f && em2 > 1.0f) {
        ik2++;
      }
      sgsum += brsig;
      // C                ratio=sgsum/ik0/40.
    }
  }
  if (icase == 2) {
    inpion++;
  }
  if (icase == 5) {
    ipipi++;
  }
  // C        write(62,*)'ik1,ik2,ik3',ik1,ik2,ik3,ratio_1,ratio_3,ratio
  // C        write(62,*)'inpion,ipipi',inpion,ipipi
  if (ranart(nseed) > (brsig / sig)) {
    // C     no anti-kaon production
    ictrl = -1;
    return;
  }
  il++;
  // C        kaons could be created now.
  if (icase == 1) {
    in++;
    // C          write(60,*)'------in,s2kaon,sig=',in,brsig,sig,lb1,lb2
    nnkaon(cmn, irun, iseed, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg);
  }
  if (icase == 2) {
    im++;
    // C          call npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
    // C     &              pcx,pcy,pcz,nchrg,ratiok)
    npik(cmn, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg,
         ratiok, iblock);
  }
  // C
  if (icase == 3) {
    im3++;
    // C          write(63,*)'im3,lb1,lb2,pkaon',im3,lb1,lb2,pkaon
    // C          write(63,*)'sig,el,sigma',sig,brel,brsgm
    // C          write(63,*)'srt,pcx,pcy,pcz,em1,em2',srt,pcx,pcy,pcz,em1,em2
    kaonn(cmn, brel, brsgm, irun, iseed, dt, nt, ictrl, i1, i2, iblock, srt,
          pcx, pcy, pcz, nchrg);
    // C         this subroutine format is diff. since three final states are
    // possible
  }
  // C
  if (icase == 4) {
    im4++;
    // C          write(64,*)'im4,sigma0,branch,sig=',im4,sigma0,brsig,sig
    // C          write(64,*)'lb1,lb2,em1,em2,pkaon=',lb1,lb2,em1,em2,pkaon
    // C
    // Csp06/07/01
    if (ranart(nseed) < brel) {
      ielstc = 1;
    } else {
      ielstc = 0;
    }
    // C          call Pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,srt,
    // C     &                   pcx,pcy,pcz,nchrg)
    pihypn(cmn, ielstc, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz,
           nchrg, iblock);
    // C
    // Csp06/07/01 end
  }
  // C        if(icase.eq.5) then
  // C          im5=im5+1
  // C          write(65,*)'---im5,s2kaon,sig=',im5,brsig,sig
  // C          call pipikaon(irun,iseed,dt,nt,ictrl,i1,i2,srt,pcx,pcy,pcz)
  // C        endif
  // Cbz3/2/99
  // C        write(101,*)lb1,lb2,lb(i1),lb(i2)
  // C        write(101,*)em1,em2,e(i1),e(i2),srt
  // Cbz3/2/99end
  // C
}

struct aknpsg_save {
  float sigma1;

  aknpsg_save() : sigma1(fem::float0) {}
};

// C
// C GQ Li parametrization (without resonance)
float aknpsg(common& cmn, float const& pkaon) {
  float return_value = fem::float0;
  FEM_CMN_SVE(aknpsg);
  // SAVE
  float& sigma1 = sve.sigma1;
  //
  // Ccross section in mb for K- + N reactions.
  // C       sigma1: x section for K- + p/n -> sigma0 + PI0
  if (pkaon <= 0.345f) {
    sigma1 = 0.624f * fem::pow(pkaon, (-1.83f));
  } else {
    sigma1 = 0.7f * fem::pow(pkaon, (-2.09f));
  }
  return_value = sigma1;
  return return_value;
}

struct artan1_save {
  float dxmt;
  float ee;
  float eta;
  int i;
  int ieta;
  int imt;
  int ityp;
  int iy;
  int j;
  float ptot;
  float px;
  float py;
  float pz;
  float xm;
  float xmt;
  float y;

  artan1_save()
      : dxmt(fem::float0),
        ee(fem::float0),
        eta(fem::float0),
        i(fem::int0),
        ieta(fem::int0),
        imt(fem::int0),
        ityp(fem::int0),
        iy(fem::int0),
        j(fem::int0),
        ptot(fem::float0),
        px(fem::float0),
        py(fem::float0),
        pz(fem::float0),
        xm(fem::float0),
        xmt(fem::float0),
        y(fem::float0) {}
};

// C
// C=======================================================================
// C
// Clin Below is the previous artana.f:
// C=======================================================================
// C
// C.....analysis subroutine before the hadronic space-time evolution
// C
void artan1(common& cmn) {
  FEM_CMN_SVE(artan1);
  common_write write(cmn);
  const int maxr = 1;
  arr_cref<int> multi1(cmn.multi1, dimension(maxr));
  const int maxstr = 150001;
  arr_cref<int, 2> ityp1(cmn.ityp1, dimension(maxstr, maxr));
  arr_cref<float, 2> px1(cmn.px1, dimension(maxstr, maxr));
  arr_cref<float, 2> py1(cmn.py1, dimension(maxstr, maxr));
  arr_cref<float, 2> pz1(cmn.pz1, dimension(maxstr, maxr));
  arr_cref<float, 2> ee1(cmn.ee1, dimension(maxstr, maxr));
  arr_cref<float, 2> xm1(cmn.xm1, dimension(maxstr, maxr));
  arr_ref<float> dy1ntb(cmn.dy1ntb, dimension(50));
  arr_ref<float> dy1ntp(cmn.dy1ntp, dimension(50));
  arr_ref<float> dy1hm(cmn.dy1hm, dimension(50));
  arr_ref<float> dy1kp(cmn.dy1kp, dimension(50));
  arr_ref<float> dy1km(cmn.dy1km, dimension(50));
  arr_ref<float> dy1k0s(cmn.dy1k0s, dimension(50));
  arr_ref<float> dy1la(cmn.dy1la, dimension(50));
  arr_ref<float> dy1lb(cmn.dy1lb, dimension(50));
  arr_ref<float> dy1phi(cmn.dy1phi, dimension(50));
  arr_ref<float> dm1pip(cmn.dm1pip, dimension(50));
  arr_ref<float> dm1pim(cmn.dm1pim, dimension(50));
  arr_ref<float> dmt1pr(cmn.dmt1pr, dimension(50));
  arr_ref<float> dmt1pb(cmn.dmt1pb, dimension(50));
  arr_ref<float> dmt1kp(cmn.dmt1kp, dimension(50));
  arr_ref<float> dm1km(cmn.dm1km, dimension(50));
  arr_ref<float> dm1k0s(cmn.dm1k0s, dimension(50));
  arr_ref<float> dmt1la(cmn.dmt1la, dimension(50));
  arr_ref<float> dmt1lb(cmn.dmt1lb, dimension(50));
  arr_ref<float> dy1msn(cmn.dy1msn, dimension(50));
  arr_ref<float> dy1pip(cmn.dy1pip, dimension(50));
  arr_ref<float> dy1pim(cmn.dy1pim, dimension(50));
  arr_ref<float> dy1pi0(cmn.dy1pi0, dimension(50));
  arr_ref<float> dy1pr(cmn.dy1pr, dimension(50));
  arr_ref<float> dy1pb(cmn.dy1pb, dimension(50));
  arr_ref<float> dy1neg(cmn.dy1neg, dimension(50));
  arr_ref<float> dy1ch(cmn.dy1ch, dimension(50));
  arr_ref<float> de1neg(cmn.de1neg, dimension(50));
  arr_ref<float> de1ch(cmn.de1ch, dimension(50));
  //
  float& dxmt = sve.dxmt;
  float& eta = sve.eta;
  int& i = sve.i;
  int& ieta = sve.ieta;
  int& imt = sve.imt;
  int& ityp = sve.ityp;
  int& iy = sve.iy;
  int& j = sve.j;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& xm = sve.xm;
  float& xmt = sve.xmt;
  float& y = sve.y;
  const float by = 0.4f;
  const float ymt1 = -1.0f;
  const float ymt2 = 1.0f;
  const float bmt = 0.05f;
  // C.....y cut for mt spectrum
  // Cbz3/17/99
  // C      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
  // Cbz3/17/99 end
  // C.....bin width for mt spectrum and y spectrum
  // Clin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
  // C      PARAMETER (BMT = 0.05, BY = 0.2)
  // Cc      SAVE /RUN/
  // Cc      SAVE /ARERC1/
  // Cbz3/17/99
  // C     &     dm1k0s(50), DMT1LA(50), DMT1LB(50)
  // Cc      SAVE /ARPRC1/
  // Cc      SAVE /ARANA1/
  // C
  // Cbz3/17/99 end
  FEM_DO_SAFE(j, 1, cmn.num) {
    FEM_DO_SAFE(i, 1, multi1(j)) {
      ityp = ityp1(i, j);
      px = px1(i, j);
      py = py1(i, j);
      pz = pz1(i, j);
      sve.ee = ee1(i, j);
      xm = xm1(i, j);
      // C     2/24/03 leptons and photons:
      if (xm < 0.01f) {
        goto statement_200;
      }
      sve.ptot = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pz));
      // Clin-9/2012 determine pseudo-rapidity more generally:
      // C            eta = 0.5*alog((Ptot+pz+1e-5)/(ptot-pz+1e-5))
      if ((fem::pow2(px) + fem::pow2(py)) > 0.f) {
        eta = asinh(pz / fem::sqrt(fem::pow2(px) + fem::pow2(py)));
      } else {
        eta = 1000000.0f * fem::sign(1.f, pz);
        // Clin-2/2013 for spectator target nucleons in LAB frame,
        // C     note that finite precision of HBOOST
        // C     would give spectator target nucleons a small but non-zero pz:
        if (fem::abs(pz) <= 1e-3f) {
          eta = 0.f;
        }
      }
      // C
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(xm));
      dxmt = xmt - xm;
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. EE) THEN
      // C               PRINT *, 'IN ARTAN1'
      // C               PRINT *, 'PARTICLE ', I, ' RUN ', J, 'PREC ERR'
      // Ccbzdbg2/16/99
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // Ccbzdbg2/16/99
      // Ccbzdbg2/15/99
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', EE
      // Ccbzdbg2/16/99
      // C               PRINT *, ' XM = ', XM
      // Ccbzdbg2/16/99end
      // C               GOTO 200
      // C            else
      // Cc            Y = 0.5 * LOG((EE + PZ) / (EE - PZ))
      // C               Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ +1e-5))
      // Cc               STOP
      // Ccbzdbg2/15/99end
      // C            END IF
      if (xmt > 0.f) {
        y = asinh(pz / xmt);
      } else {
        write(6, star), " IN ARTAN1 mt=0";
        y = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      // C.....rapidity cut for the rapidity distribution
      if (fem::abs(y) >= 10.0f) {
        goto statement_100;
      }
      // Clin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
      // C            IY = 1 + int(ABS(Y) / BY)
      // C            Ieta = 1 + int(ABS(eta) / BY)
      if (fem::abs(eta) >= 10.0f) {
        goto statement_100;
      }
      iy = 1 + fem::fint((y + 10.f) / by);
      ieta = 1 + fem::fint((eta + 10.f) / by);
      // C
      if (ityp < -1000) {
        dy1ntb(iy) = dy1ntb(iy) - 1.0f;
      }
      if (ityp > 1000) {
        dy1ntb(iy) += 1.0f;
      }
      if (ityp == -2212) {
        dy1ntp(iy) = dy1ntp(iy) - 1.0f;
      }
      if (ityp == 2212) {
        dy1ntp(iy) += 1.0f;
      }
      // C            IF (ITYP .EQ. -211 .OR. ITYP .EQ. -321 .OR.
      // C     &         ITYP .EQ. -2212) THEN
      if (ityp == -2112) {
        dy1hm(iy) += 1.0f;
      }
      // C
      if (luchge(ityp) != 0) {
        dy1ch(iy) += 1.0f;
        de1ch(ieta) += 1.0f;
        if (luchge(ityp) < 0) {
          dy1neg(iy) += 1.0f;
          de1neg(ieta) += 1.0f;
        }
      }
      // C
      // Cbz3/17/99
      if ((ityp >= 100 && ityp < 1000) || (ityp > -1000 && ityp <= -100)) {
        dy1msn(iy) += 1.0f;
      }
      if (ityp == 211) {
        dy1pip(iy) += 1.0f;
      }
      if (ityp == -211) {
        dy1pim(iy) += 1.0f;
      }
      if (ityp == 111) {
        dy1pi0(iy) += 1.0f;
      }
      if (ityp == 2212) {
        dy1pr(iy) += 1.0f;
      }
      if (ityp == -2212) {
        dy1pb(iy) += 1.0f;
      }
      // Cbz3/17/99 end
      if (ityp == 321) {
        dy1kp(iy) += 1.0f;
      }
      if (ityp == -321) {
        dy1km(iy) += 1.0f;
      }
      // Clin-4/24/03 evaluate K0L instead of K0S, since sometimes we may decay
      // K0S: C            IF (ITYP .EQ. 310) THEN
      if (ityp == 130) {
        dy1k0s(iy) += 1.0f;
      }
      if (ityp == 3122) {
        dy1la(iy) += 1.0f;
      }
      if (ityp == -3122) {
        dy1lb(iy) += 1.0f;
      }
      if (ityp == 333) {
        dy1phi(iy) += 1.0f;
      }
    // C
    // C.....insert rapidity cut for mt spectrum here
    statement_100:
      if (y < ymt1 || y > ymt2) {
        goto statement_200;
      }
      if (dxmt >= 50.0f * bmt || dxmt == 0) {
        goto statement_200;
      }
      imt = 1 + fem::fint(dxmt / bmt);
      if (ityp == 211) {
        dm1pip(imt) += 1.0f / xmt;
      }
      if (ityp == -211) {
        dm1pim(imt) += 1.0f / xmt;
      }
      if (ityp == 2212) {
        dmt1pr(imt) += 1.0f / xmt;
      }
      if (ityp == -2212) {
        dmt1pb(imt) += 1.0f / xmt;
      }
      if (ityp == 321) {
        dmt1kp(imt) += 1.0f / xmt;
      }
      if (ityp == -321) {
        dm1km(imt) += 1.0f / xmt;
      }
      // Clin-4/24/03:
      // C            IF (ITYP .EQ. 310) THEN
      if (ityp == 130) {
        dm1k0s(imt) += 1.0f / xmt;
      }
      if (ityp == 3122) {
        dmt1la(imt) += 1.0f / xmt;
      }
      if (ityp == -3122) {
        dmt1lb(imt) += 1.0f / xmt;
      }
    // C
    statement_200:;
    }
  }
  // C
}

struct artan2_save {
  float dxmt;
  float ee;
  float eta;
  int i;
  int ieta;
  int imt;
  int ityp;
  int iy;
  int j;
  float ptot;
  float px;
  float py;
  float pz;
  float xm;
  float xmt;
  float y;

  artan2_save()
      : dxmt(fem::float0),
        ee(fem::float0),
        eta(fem::float0),
        i(fem::int0),
        ieta(fem::int0),
        imt(fem::int0),
        ityp(fem::int0),
        iy(fem::int0),
        j(fem::int0),
        ptot(fem::float0),
        px(fem::float0),
        py(fem::float0),
        pz(fem::float0),
        xm(fem::float0),
        xmt(fem::float0),
        y(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine after the hadronic space-time evolution
// C
void artan2(common& cmn) {
  FEM_CMN_SVE(artan2);
  common_write write(cmn);
  const int maxr = 1;
  arr_cref<int> multi1(cmn.multi1, dimension(maxr));
  const int maxstr = 150001;
  arr_cref<int, 2> ityp1(cmn.ityp1, dimension(maxstr, maxr));
  arr_cref<float, 2> px1(cmn.px1, dimension(maxstr, maxr));
  arr_cref<float, 2> py1(cmn.py1, dimension(maxstr, maxr));
  arr_cref<float, 2> pz1(cmn.pz1, dimension(maxstr, maxr));
  arr_cref<float, 2> ee1(cmn.ee1, dimension(maxstr, maxr));
  arr_cref<float, 2> xm1(cmn.xm1, dimension(maxstr, maxr));
  arr_ref<float> dy2ntb(cmn.dy2ntb, dimension(50));
  arr_ref<float> dy2ntp(cmn.dy2ntp, dimension(50));
  arr_ref<float> dy2hm(cmn.dy2hm, dimension(50));
  arr_ref<float> dy2kp(cmn.dy2kp, dimension(50));
  arr_ref<float> dy2km(cmn.dy2km, dimension(50));
  arr_ref<float> dy2k0s(cmn.dy2k0s, dimension(50));
  arr_ref<float> dy2la(cmn.dy2la, dimension(50));
  arr_ref<float> dy2lb(cmn.dy2lb, dimension(50));
  arr_ref<float> dy2phi(cmn.dy2phi, dimension(50));
  arr_ref<float> dm2pip(cmn.dm2pip, dimension(50));
  arr_ref<float> dm2pim(cmn.dm2pim, dimension(50));
  arr_ref<float> dmt2pr(cmn.dmt2pr, dimension(50));
  arr_ref<float> dmt2pb(cmn.dmt2pb, dimension(50));
  arr_ref<float> dmt2kp(cmn.dmt2kp, dimension(50));
  arr_ref<float> dm2km(cmn.dm2km, dimension(50));
  arr_ref<float> dm2k0s(cmn.dm2k0s, dimension(50));
  arr_ref<float> dmt2la(cmn.dmt2la, dimension(50));
  arr_ref<float> dmt2lb(cmn.dmt2lb, dimension(50));
  arr_ref<float> dy2msn(cmn.dy2msn, dimension(50));
  arr_ref<float> dy2pip(cmn.dy2pip, dimension(50));
  arr_ref<float> dy2pim(cmn.dy2pim, dimension(50));
  arr_ref<float> dy2pi0(cmn.dy2pi0, dimension(50));
  arr_ref<float> dy2pr(cmn.dy2pr, dimension(50));
  arr_ref<float> dy2pb(cmn.dy2pb, dimension(50));
  arr_ref<float> dy2neg(cmn.dy2neg, dimension(50));
  arr_ref<float> dy2ch(cmn.dy2ch, dimension(50));
  arr_ref<float> de2neg(cmn.de2neg, dimension(50));
  arr_ref<float> de2ch(cmn.de2ch, dimension(50));
  //
  float& dxmt = sve.dxmt;
  float& eta = sve.eta;
  int& i = sve.i;
  int& ieta = sve.ieta;
  int& imt = sve.imt;
  int& ityp = sve.ityp;
  int& iy = sve.iy;
  int& j = sve.j;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& xm = sve.xm;
  float& xmt = sve.xmt;
  float& y = sve.y;
  const float by = 0.4f;
  const float ymt1 = -1.0f;
  const float ymt2 = 1.0f;
  const float bmt = 0.05f;
  // C
  // C.....y cut for mt spectrum
  // Cbz3/17/99
  // C      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
  // Cbz3/17/99 end
  // C.....bin width for mt spectrum and y spectrum
  // C      PARAMETER (BMT = 0.05, BY = 0.2)
  // Cc      SAVE /RUN/
  // Cc      SAVE /ARERC1/
  // Cbz3/17/99
  // C     &     dm2k0s(50), DMT2LA(50), DMT2LB(50)
  // Cc      SAVE /ARPRC1/
  // Cbz3/17/99 end
  // Cc      SAVE /ARANA2/
  // C
  FEM_DO_SAFE(j, 1, cmn.num) {
    FEM_DO_SAFE(i, 1, multi1(j)) {
      ityp = ityp1(i, j);
      px = px1(i, j);
      py = py1(i, j);
      pz = pz1(i, j);
      sve.ee = ee1(i, j);
      xm = xm1(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(xm));
      // C     2/24/03 leptons and photons:
      if (xm < 0.01f) {
        goto statement_200;
      }
      dxmt = xmt - xm;
      sve.ptot = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pz));
      // Clin-9/2012 determine pseudo-rapidity more generally:
      // C            eta = 0.5*alog((Ptot+pz+1e-5)/(ptot-pz+1e-5))
      if ((fem::pow2(px) + fem::pow2(py)) > 0.f) {
        eta = asinh(pz / fem::sqrt(fem::pow2(px) + fem::pow2(py)));
      } else {
        eta = 1000000.0f * fem::sign(1.f, pz);
        // Clin-2/2013 for spectator target nucleons in LAB frame:
        if (fem::abs(pz) <= 1e-3f) {
          eta = 0.f;
        }
      }
      // C
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. EE) THEN
      // C               PRINT *, 'IN ARTAN2'
      // C               PRINT *, 'PARTICLE ', I, ' RUN ', J, 'PREC ERR'
      // Ccbzdbg2/16/99
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // Ccbzdbg2/16/99
      // Ccbzdbg2/15/99
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', EE
      // Ccbzdbg2/16/99
      // C               PRINT *, ' XM = ', XM
      // Ccbzdbg2/16/99end
      // C               GOTO 200
      // Cc               STOP
      // Ccbzdbg2/15/99end
      // C            END IF
      // C            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
      if (xmt > 0.f) {
        y = asinh(pz / xmt);
      } else {
        write(6, star), " IN ARTAN2 mt=0";
        y = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      // C.....rapidity cut for the rapidity distribution
      if (fem::abs(y) >= 10.0f) {
        goto statement_100;
      }
      // C            IY = 1 + int(ABS(Y) / BY)
      // C            Ieta = 1 + int(ABS(eta) / BY)
      if (fem::abs(eta) >= 10.0f) {
        goto statement_100;
      }
      iy = 1 + fem::fint((y + 10.f) / by);
      ieta = 1 + fem::fint((eta + 10.f) / by);
      // C
      if (ityp < -1000) {
        dy2ntb(iy) = dy2ntb(iy) - 1.0f;
      }
      if (ityp > 1000) {
        dy2ntb(iy) += 1.0f;
      }
      if (ityp == -2212) {
        dy2ntp(iy) = dy2ntp(iy) - 1.0f;
      }
      if (ityp == 2212) {
        dy2ntp(iy) += 1.0f;
      }
      if (ityp == -2112) {
        dy2hm(iy) += 1.0f;
      }
      // C
      if (luchge(ityp) != 0) {
        dy2ch(iy) += 1.0f;
        de2ch(ieta) += 1.0f;
        if (luchge(ityp) < 0) {
          dy2neg(iy) += 1.0f;
          de2neg(ieta) += 1.0f;
        }
      }
      // C
      // Cbz3/17/99
      if ((ityp >= 100 && ityp < 1000) || (ityp > -1000 && ityp <= -100)) {
        dy2msn(iy) += 1.0f;
      }
      if (ityp == 211) {
        dy2pip(iy) += 1.0f;
      }
      if (ityp == -211) {
        dy2pim(iy) += 1.0f;
      }
      if (ityp == 111) {
        dy2pi0(iy) += 1.0f;
      }
      if (ityp == 2212) {
        dy2pr(iy) += 1.0f;
      }
      if (ityp == -2212) {
        dy2pb(iy) += 1.0f;
      }
      // Cbz3/17/99 end
      if (ityp == 321) {
        dy2kp(iy) += 1.0f;
      }
      if (ityp == -321) {
        dy2km(iy) += 1.0f;
      }
      // Clin-4/24/03:
      // C            IF (ITYP .EQ. 310) THEN
      if (ityp == 130) {
        dy2k0s(iy) += 1.0f;
      }
      if (ityp == 3122) {
        dy2la(iy) += 1.0f;
      }
      if (ityp == -3122) {
        dy2lb(iy) += 1.0f;
      }
      if (ityp == 333) {
        dy2phi(iy) += 1.0f;
      }
    // C
    // C.....insert rapidity cut for mt spectrum here
    statement_100:
      if (y < ymt1 || y > ymt2) {
        goto statement_200;
      }
      if (dxmt >= 50.0f * bmt || dxmt == 0) {
        goto statement_200;
      }
      imt = 1 + fem::fint(dxmt / bmt);
      if (ityp == 211) {
        dm2pip(imt) += 1.0f / xmt;
      }
      if (ityp == -211) {
        dm2pim(imt) += 1.0f / xmt;
      }
      if (ityp == 2212) {
        dmt2pr(imt) += 1.0f / xmt;
      }
      if (ityp == -2212) {
        dmt2pb(imt) += 1.0f / xmt;
      }
      if (ityp == 321) {
        dmt2kp(imt) += 1.0f / xmt;
      }
      if (ityp == -321) {
        dm2km(imt) += 1.0f / xmt;
      }
      // Clin-4/24/03:
      // C            IF (ITYP .EQ. 310) THEN
      if (ityp == 130) {
        dm2k0s(imt) += 1.0f / xmt;
      }
      if (ityp == 3122) {
        dmt2la(imt) += 1.0f / xmt;
      }
      if (ityp == -3122) {
        dmt2lb(imt) += 1.0f / xmt;
      }
    // C
    statement_200:;
    }
  }
  // C
}

struct artout_save {
  int i;
  float scale1;
  float scale2;
  float ymid;

  artout_save()
      : i(fem::int0),
        scale1(fem::float0),
        scale2(fem::float0),
        ymid(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....output analysis results at the end of the simulation
// C
void artout(common& cmn, int const& nevnt) {
  FEM_CMN_SVE(artout);
  common_write write(cmn);
  // COMMON run
  int& num = cmn.num;
  // COMMON arana1
  arr_cref<float> dy1ntb(cmn.dy1ntb, dimension(50));
  arr_cref<float> dy1ntp(cmn.dy1ntp, dimension(50));
  arr_cref<float> dy1hm(cmn.dy1hm, dimension(50));
  arr_cref<float> dy1kp(cmn.dy1kp, dimension(50));
  arr_cref<float> dy1km(cmn.dy1km, dimension(50));
  arr_cref<float> dy1k0s(cmn.dy1k0s, dimension(50));
  arr_cref<float> dy1la(cmn.dy1la, dimension(50));
  arr_cref<float> dy1lb(cmn.dy1lb, dimension(50));
  arr_cref<float> dy1phi(cmn.dy1phi, dimension(50));
  arr_cref<float> dm1pip(cmn.dm1pip, dimension(50));
  arr_cref<float> dm1pim(cmn.dm1pim, dimension(50));
  arr_cref<float> dmt1pr(cmn.dmt1pr, dimension(50));
  arr_cref<float> dmt1pb(cmn.dmt1pb, dimension(50));
  arr_cref<float> dmt1kp(cmn.dmt1kp, dimension(50));
  arr_cref<float> dm1km(cmn.dm1km, dimension(50));
  arr_cref<float> dm1k0s(cmn.dm1k0s, dimension(50));
  arr_cref<float> dmt1la(cmn.dmt1la, dimension(50));
  arr_cref<float> dmt1lb(cmn.dmt1lb, dimension(50));
  arr_cref<float> dy1msn(cmn.dy1msn, dimension(50));
  arr_cref<float> dy1pip(cmn.dy1pip, dimension(50));
  arr_cref<float> dy1pim(cmn.dy1pim, dimension(50));
  arr_cref<float> dy1pi0(cmn.dy1pi0, dimension(50));
  arr_cref<float> dy1pr(cmn.dy1pr, dimension(50));
  arr_cref<float> dy1pb(cmn.dy1pb, dimension(50));
  arr_cref<float> dy1neg(cmn.dy1neg, dimension(50));
  arr_cref<float> dy1ch(cmn.dy1ch, dimension(50));
  arr_cref<float> de1neg(cmn.de1neg, dimension(50));
  arr_cref<float> de1ch(cmn.de1ch, dimension(50));
  // COMMON arana2
  arr_cref<float> dy2ntb(cmn.dy2ntb, dimension(50));
  arr_cref<float> dy2ntp(cmn.dy2ntp, dimension(50));
  arr_cref<float> dy2hm(cmn.dy2hm, dimension(50));
  arr_cref<float> dy2kp(cmn.dy2kp, dimension(50));
  arr_cref<float> dy2km(cmn.dy2km, dimension(50));
  arr_cref<float> dy2k0s(cmn.dy2k0s, dimension(50));
  arr_cref<float> dy2la(cmn.dy2la, dimension(50));
  arr_cref<float> dy2lb(cmn.dy2lb, dimension(50));
  arr_cref<float> dy2phi(cmn.dy2phi, dimension(50));
  arr_cref<float> dm2pip(cmn.dm2pip, dimension(50));
  arr_cref<float> dm2pim(cmn.dm2pim, dimension(50));
  arr_cref<float> dmt2pr(cmn.dmt2pr, dimension(50));
  arr_cref<float> dmt2pb(cmn.dmt2pb, dimension(50));
  arr_cref<float> dmt2kp(cmn.dmt2kp, dimension(50));
  arr_cref<float> dm2km(cmn.dm2km, dimension(50));
  arr_cref<float> dm2k0s(cmn.dm2k0s, dimension(50));
  arr_cref<float> dmt2la(cmn.dmt2la, dimension(50));
  arr_cref<float> dmt2lb(cmn.dmt2lb, dimension(50));
  arr_cref<float> dy2msn(cmn.dy2msn, dimension(50));
  arr_cref<float> dy2pip(cmn.dy2pip, dimension(50));
  arr_cref<float> dy2pim(cmn.dy2pim, dimension(50));
  arr_cref<float> dy2pi0(cmn.dy2pi0, dimension(50));
  arr_cref<float> dy2pr(cmn.dy2pr, dimension(50));
  arr_cref<float> dy2pb(cmn.dy2pb, dimension(50));
  arr_cref<float> dy2neg(cmn.dy2neg, dimension(50));
  arr_cref<float> dy2ch(cmn.dy2ch, dimension(50));
  arr_cref<float> de2neg(cmn.de2neg, dimension(50));
  arr_cref<float> de2ch(cmn.de2ch, dimension(50));
  //
  // SAVE
  int& i = sve.i;
  float& scale1 = sve.scale1;
  float& scale2 = sve.scale2;
  float& ymid = sve.ymid;
  //
  static const char* format_333 = "(2(f12.5,1x))";
  // C
  // C.....y cut for mt spectrum
  // Cbz3/17/99
  // C      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
  // Cbz3/17/99 end
  // C.....bin width for mt spectrum and y spectrum
  // C      PARAMETER (BMT = 0.05, BY = 0.2)
  // Cc      SAVE /RUN/
  // Cbz3/17/99
  // C     &     dm1k0s(50), DMT1LA(50), DMT1LB(50)
  // Cc      SAVE /ARPRC1/
  // Cbz3/17/99 end
  // Cc      SAVE /ARANA1/
  // Cbz3/17/99
  // C     &     dm2k0s(50), DMT2LA(50), DMT2LB(50)
  // Cc      SAVE /ARANA2/
  // Cbz3/17/99 end
  cmn.io.open(30, "ana/dndy_netb.dat").status("UNKNOWN");
  cmn.io.open(31, "ana/dndy_netp.dat").status("UNKNOWN");
  cmn.io.open(32, "ana/dndy_nb.dat").status("UNKNOWN");
  cmn.io.open(33, "ana/dndy_neg.dat").status("UNKNOWN");
  cmn.io.open(34, "ana/dndy_ch.dat").status("UNKNOWN");
  cmn.io.open(35, "ana/dnde_neg.dat").status("UNKNOWN");
  cmn.io.open(36, "ana/dnde_ch.dat").status("UNKNOWN");
  cmn.io.open(37, "ana/dndy_kp.dat").status("UNKNOWN");
  cmn.io.open(38, "ana/dndy_km.dat").status("UNKNOWN");
  // Clin-4/24/03
  // C      OPEN (39, FILE = 'ana/dndy_k0s.dat', STATUS = 'UNKNOWN')
  cmn.io.open(39, "ana/dndy_k0l.dat").status("UNKNOWN");
  cmn.io.open(40, "ana/dndy_la.dat").status("UNKNOWN");
  cmn.io.open(41, "ana/dndy_lb.dat").status("UNKNOWN");
  cmn.io.open(42, "ana/dndy_phi.dat").status("UNKNOWN");
  // Cbz3/17/99
  cmn.io.open(43, "ana/dndy_meson.dat").status("UNKNOWN");
  cmn.io.open(44, "ana/dndy_pip.dat").status("UNKNOWN");
  cmn.io.open(45, "ana/dndy_pim.dat").status("UNKNOWN");
  cmn.io.open(46, "ana/dndy_pi0.dat").status("UNKNOWN");
  cmn.io.open(47, "ana/dndy_pr.dat").status("UNKNOWN");
  cmn.io.open(48, "ana/dndy_pb.dat").status("UNKNOWN");
  // Cbz3/17/99 end
  // C
  cmn.io.open(50, "ana/dndmtdy_pip.dat").status("UNKNOWN");
  cmn.io.open(51, "ana/dndmtdy_0_1_pim.dat").status("UNKNOWN");
  cmn.io.open(52, "ana/dndmtdy_pr.dat").status("UNKNOWN");
  cmn.io.open(53, "ana/dndmtdy_pb.dat").status("UNKNOWN");
  cmn.io.open(54, "ana/dndmtdy_kp.dat").status("UNKNOWN");
  cmn.io.open(55, "ana/dndmtdy_0_5_km.dat").status("UNKNOWN");
  cmn.io.open(56, "ana/dndmtdy_k0s.dat").status("UNKNOWN");
  cmn.io.open(57, "ana/dndmtdy_la.dat").status("UNKNOWN");
  cmn.io.open(58, "ana/dndmtdy_lb.dat").status("UNKNOWN");
  // Clin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
  // C      SCALE1 = 1. / REAL(NEVNT * NUM) / BY / 2.0
  const float by = 0.4f;
  scale1 = 1.f / fem::real(nevnt * num) / by;
  const float bmt = 0.05f;
  const float ymt2 = 1.0f;
  const float ymt1 = -1.0f;
  scale2 = 1.f / fem::real(nevnt * num) / bmt / (ymt2 - ymt1);
  // C
  FEM_DO_SAFE(i, 1, 50) {
    ymid = -10.f + by * (i - 0.5f);
    write(30, format_333), ymid, scale1* dy1ntb(i);
    write(31, format_333), ymid, scale1* dy1ntp(i);
    write(32, format_333), ymid, scale1* dy1hm(i);
    write(37, format_333), ymid, scale1* dy1kp(i);
    write(38, format_333), ymid, scale1* dy1km(i);
    write(39, format_333), ymid, scale1* dy1k0s(i);
    write(40, format_333), ymid, scale1* dy1la(i);
    write(41, format_333), ymid, scale1* dy1lb(i);
    write(42, format_333), ymid, scale1* dy1phi(i);
    write(33, format_333), ymid, scale1* dy1neg(i);
    write(34, format_333), ymid, scale1* dy1ch(i);
    write(35, format_333), ymid, scale1* de1neg(i);
    write(36, format_333), ymid, scale1* de1ch(i);
    write(43, format_333), ymid, scale1* dy1msn(i);
    write(44, format_333), ymid, scale1* dy1pip(i);
    write(45, format_333), ymid, scale1* dy1pim(i);
    write(46, format_333), ymid, scale1* dy1pi0(i);
    write(47, format_333), ymid, scale1* dy1pr(i);
    write(48, format_333), ymid, scale1* dy1pb(i);
    // C
    if (dm1pip(i) != 0.0f) {
      write(50, format_333), bmt*(i - 0.5f), scale2* dm1pip(i);
    }
    if (dm1pim(i) != 0.0f) {
      write(51, format_333), bmt*(i - 0.5f), scale2 * 0.1f * dm1pim(i);
    }
    if (dmt1pr(i) != 0.0f) {
      write(52, format_333), bmt*(i - 0.5f), scale2* dmt1pr(i);
    }
    if (dmt1pb(i) != 0.0f) {
      write(53, format_333), bmt*(i - 0.5f), scale2* dmt1pb(i);
    }
    if (dmt1kp(i) != 0.0f) {
      write(54, format_333), bmt*(i - 0.5f), scale2* dmt1kp(i);
    }
    if (dm1km(i) != 0.0f) {
      write(55, format_333), bmt*(i - 0.5f), scale2 * 0.5f * dm1km(i);
    }
    if (dm1k0s(i) != 0.0f) {
      write(56, format_333), bmt*(i - 0.5f), scale2* dm1k0s(i);
    }
    if (dmt1la(i) != 0.0f) {
      write(57, format_333), bmt*(i - 0.5f), scale2* dmt1la(i);
    }
    if (dmt1lb(i) != 0.0f) {
      write(58, format_333), bmt*(i - 0.5f), scale2* dmt1lb(i);
    }
  }
  // C
  FEM_DO_SAFE(i, 30, 48) { write(i, star), "after hadron evolution"; }
  FEM_DO_SAFE(i, 50, 58) { write(i, star), "after hadron evolution"; }
  // C
  FEM_DO_SAFE(i, 1, 50) {
    ymid = -10.f + by * (i - 0.5f);
    write(30, format_333), ymid, scale1* dy2ntb(i);
    write(31, format_333), ymid, scale1* dy2ntp(i);
    write(32, format_333), ymid, scale1* dy2hm(i);
    write(37, format_333), ymid, scale1* dy2kp(i);
    write(38, format_333), ymid, scale1* dy2km(i);
    write(39, format_333), ymid, scale1* dy2k0s(i);
    write(40, format_333), ymid, scale1* dy2la(i);
    write(41, format_333), ymid, scale1* dy2lb(i);
    write(42, format_333), ymid, scale1* dy2phi(i);
    write(33, format_333), ymid, scale1* dy2neg(i);
    write(34, format_333), ymid, scale1* dy2ch(i);
    write(35, format_333), ymid, scale1* de2neg(i);
    write(36, format_333), ymid, scale1* de2ch(i);
    write(43, format_333), ymid, scale1* dy2msn(i);
    write(44, format_333), ymid, scale1* dy2pip(i);
    write(45, format_333), ymid, scale1* dy2pim(i);
    write(46, format_333), ymid, scale1* dy2pi0(i);
    write(47, format_333), ymid, scale1* dy2pr(i);
    write(48, format_333), ymid, scale1* dy2pb(i);
    // C
    if (dm2pip(i) != 0.0f) {
      write(50, format_333), bmt*(i - 0.5f), scale2* dm2pip(i);
    }
    if (dm2pim(i) != 0.0f) {
      write(51, format_333), bmt*(i - 0.5f), scale2 * 0.1f * dm2pim(i);
    }
    if (dmt2pr(i) != 0.0f) {
      write(52, format_333), bmt*(i - 0.5f), scale2* dmt2pr(i);
    }
    if (dmt2pb(i) != 0.0f) {
      write(53, format_333), bmt*(i - 0.5f), scale2* dmt2pb(i);
    }
    if (dmt2kp(i) != 0.0f) {
      write(54, format_333), bmt*(i - 0.5f), scale2* dmt2kp(i);
    }
    if (dm2km(i) != 0.0f) {
      write(55, format_333), bmt*(i - 0.5f), scale2 * 0.5f * dm2km(i);
    }
    if (dm2k0s(i) != 0.0f) {
      write(56, format_333), bmt*(i - 0.5f), scale2* dm2k0s(i);
    }
    if (dmt2la(i) != 0.0f) {
      write(57, format_333), bmt*(i - 0.5f), scale2* dmt2la(i);
    }
    if (dmt2lb(i) != 0.0f) {
      write(58, format_333), bmt*(i - 0.5f), scale2* dmt2lb(i);
    }
  }
  // C
}

struct hjana1_save {
  arr<float> deyg1;
  arr<float> deyg1c;
  arr<float> deyp1;
  arr<float> dmyg1;
  arr<float> dmyg1c;
  arr<float> dmyp1;
  arr<float> dnrin1;
  arr<float> dnrpj1;
  arr<float> dnrtg1;
  arr<float> dnrtt1;
  float dxmt;
  arr<float> dyg1;
  arr<float> dyg1c;
  arr<float> dyp1;
  int i;
  int imt;
  int ir;
  int isevt;
  int isrun;
  int ityp;
  int iw;
  int iy;
  int j;
  int nisg;
  int nisgs;
  int nsubg;
  int nsubgs;
  int nsubp;
  int nsubps;
  float pe;
  float pm;
  float px;
  float py;
  float pz;
  float rap;
  arr<float> seyg1;
  arr<float> seyg1c;
  arr<float> seyp1;
  arr<float> smyg1;
  arr<float> smyg1c;
  arr<float> smyp1;
  arr<float> snrin1;
  arr<float> snrpj1;
  arr<float> snrtg1;
  arr<float> snrtt1;
  arr<float> snyg1;
  arr<float> snyg1c;
  arr<float> snyp1;
  float xmt;
  float y1;
  float y2;
  float yr;

  hjana1_save()
      : deyg1(dimension(50), fem::fill0),
        deyg1c(dimension(50), fem::fill0),
        deyp1(dimension(50), fem::fill0),
        dmyg1(dimension(200), fem::fill0),
        dmyg1c(dimension(50), fem::fill0),
        dmyp1(dimension(200), fem::fill0),
        dnrin1(dimension(50), fem::fill0),
        dnrpj1(dimension(50), fem::fill0),
        dnrtg1(dimension(50), fem::fill0),
        dnrtt1(dimension(50), fem::fill0),
        dxmt(fem::float0),
        dyg1(dimension(50), fem::fill0),
        dyg1c(dimension(50), fem::fill0),
        dyp1(dimension(50), fem::fill0),
        i(fem::int0),
        imt(fem::int0),
        ir(fem::int0),
        isevt(fem::int0),
        isrun(fem::int0),
        ityp(fem::int0),
        iw(fem::int0),
        iy(fem::int0),
        j(fem::int0),
        nisg(fem::int0),
        nisgs(fem::int0),
        nsubg(fem::int0),
        nsubgs(fem::int0),
        nsubp(fem::int0),
        nsubps(fem::int0),
        pe(fem::float0),
        pm(fem::float0),
        px(fem::float0),
        py(fem::float0),
        pz(fem::float0),
        rap(fem::float0),
        seyg1(dimension(50), fem::fill0),
        seyg1c(dimension(50), fem::fill0),
        seyp1(dimension(50), fem::fill0),
        smyg1(dimension(200), fem::fill0),
        smyg1c(dimension(50), fem::fill0),
        smyp1(dimension(200), fem::fill0),
        snrin1(dimension(50), fem::fill0),
        snrpj1(dimension(50), fem::fill0),
        snrtg1(dimension(50), fem::fill0),
        snrtt1(dimension(50), fem::fill0),
        snyg1(dimension(50), fem::fill0),
        snyg1c(dimension(50), fem::fill0),
        snyp1(dimension(50), fem::fill0),
        xmt(fem::float0),
        y1(fem::float0),
        y2(fem::float0),
        yr(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine in HIJING before parton cascade evolution
void hjana1(common& cmn) {
  FEM_CMN_SVE(hjana1);
  common_write write(cmn);
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
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
  int& nsg = cmn.nsg;
  const int maxstr = 150001;
  arr_cref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_cref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_cref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_cref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_cref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_cref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_cref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_cref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  const int maxptn = 400001;
  arr_cref<double> px0(cmn.px0, dimension(maxptn));
  arr_cref<double> py0(cmn.py0, dimension(maxptn));
  arr_cref<double> pz0(cmn.pz0, dimension(maxptn));
  arr_cref<double> e0(cmn.e0, dimension(maxptn));
  arr_cref<double> xmass0(cmn.xmass0, dimension(maxptn));
  arr_cref<int> ityp0(cmn.ityp0, dimension(maxptn));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  //
  arr_ref<float> deyg1(sve.deyg1, dimension(50));
  arr_ref<float> deyg1c(sve.deyg1c, dimension(50));
  arr_ref<float> deyp1(sve.deyp1, dimension(50));
  arr_ref<float> dmyg1(sve.dmyg1, dimension(200));
  arr_ref<float> dmyg1c(sve.dmyg1c, dimension(50));
  arr_ref<float> dmyp1(sve.dmyp1, dimension(200));
  arr_ref<float> dnrin1(sve.dnrin1, dimension(50));
  arr_ref<float> dnrpj1(sve.dnrpj1, dimension(50));
  arr_ref<float> dnrtg1(sve.dnrtg1, dimension(50));
  arr_ref<float> dnrtt1(sve.dnrtt1, dimension(50));
  float& dxmt = sve.dxmt;
  arr_ref<float> dyg1(sve.dyg1, dimension(50));
  arr_ref<float> dyg1c(sve.dyg1c, dimension(50));
  arr_ref<float> dyp1(sve.dyp1, dimension(50));
  int& i = sve.i;
  int& imt = sve.imt;
  int& ir = sve.ir;
  int& isevt = sve.isevt;
  int& isrun = sve.isrun;
  int& ityp = sve.ityp;
  int& iw = sve.iw;
  int& iy = sve.iy;
  int& j = sve.j;
  int& nisg = sve.nisg;
  int& nisgs = sve.nisgs;
  int& nsubg = sve.nsubg;
  int& nsubgs = sve.nsubgs;
  int& nsubp = sve.nsubp;
  int& nsubps = sve.nsubps;
  float& pe = sve.pe;
  float& pm = sve.pm;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& rap = sve.rap;
  arr_ref<float> seyg1(sve.seyg1, dimension(50));
  arr_ref<float> seyg1c(sve.seyg1c, dimension(50));
  arr_ref<float> seyp1(sve.seyp1, dimension(50));
  arr_ref<float> smyg1(sve.smyg1, dimension(200));
  arr_ref<float> smyg1c(sve.smyg1c, dimension(50));
  arr_ref<float> smyp1(sve.smyp1, dimension(200));
  arr_ref<float> snrin1(sve.snrin1, dimension(50));
  arr_ref<float> snrpj1(sve.snrpj1, dimension(50));
  arr_ref<float> snrtg1(sve.snrtg1, dimension(50));
  arr_ref<float> snrtt1(sve.snrtt1, dimension(50));
  arr_ref<float> snyg1(sve.snyg1, dimension(50));
  arr_ref<float> snyg1c(sve.snyg1c, dimension(50));
  arr_ref<float> snyp1(sve.snyp1, dimension(50));
  float& xmt = sve.xmt;
  float& y1 = sve.y1;
  float& y2 = sve.y2;
  float& yr = sve.yr;
  if (is_called_first_time) {
    iw = 0;
  }
  const float dy = 0.2f;
  const float dmt = 0.05f;
  const float ymax = 1.0f;
  const float ymin = -1.0f;
  const float dr = 0.2f;
  // C
  // Cc      SAVE /PARA1/
  // Cc      SAVE /HPARNT/
  // Cc      SAVE /hjcrdn/
  // Cc      SAVE /HJJET1/
  // Cc      SAVE /HJJET2/
  // Cc      SAVE /prec1/
  // Cc      SAVE /AREVT/
  // Cc      SAVE /AROUT/
  // C
  if (isevt == iaevt && isrun == iarun) {
    FEM_DO_SAFE(i, 1, 200) {
      dmyp1(i) = smyp1(i);
      dmyg1(i) = smyg1(i);
    }
    // C
    FEM_DO_SAFE(i, 1, 50) {
      dyp1(i) = snyp1(i);
      deyp1(i) = seyp1(i);
      dyg1(i) = snyg1(i);
      deyg1(i) = seyg1(i);
      dnrpj1(i) = snrpj1(i);
      dnrtg1(i) = snrtg1(i);
      dnrin1(i) = snrin1(i);
      dnrtt1(i) = snrtt1(i);
      dyg1c(i) = snyg1c(i);
      dmyg1c(i) = smyg1c(i);
      deyg1c(i) = seyg1c(i);
    }
    nsubp = nsubps;
    nsubg = nsubgs;
    nisg = nisgs;
  } else {
    FEM_DO_SAFE(i, 1, 200) {
      smyp1(i) = dmyp1(i);
      smyg1(i) = dmyg1(i);
    }
    // C
    FEM_DO_SAFE(i, 1, 50) {
      snyp1(i) = dyp1(i);
      seyp1(i) = deyp1(i);
      snyg1(i) = dyg1(i);
      seyg1(i) = deyg1(i);
      snrpj1(i) = dnrpj1(i);
      snrtg1(i) = dnrtg1(i);
      snrin1(i) = dnrin1(i);
      snrtt1(i) = dnrtt1(i);
      snyg1c(i) = dyg1c(i);
      smyg1c(i) = dmyg1c(i);
      seyg1c(i) = deyg1c(i);
    }
    nsubps = nsubp;
    nsubgs = nsubg;
    nisgs = nisg;
    isevt = iaevt;
    isrun = iarun;
    iw++;
  }
  // C.....analysis
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      ityp = kfpj(i, j);
      px = pjpx(i, j);
      py = pjpy(i, j);
      pz = pjpz(i, j);
      pe = pjpe(i, j);
      pm = pjpm(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
      dxmt = xmt - pm;
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. PE) THEN
      // C               PRINT *, ' IN HJANA1, PROJ STR ', I, ' PART ', J
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', PE
      // C               PRINT *, ' XM = ', PM
      // C               GOTO 200
      // C            END IF
      // C            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      if (xmt > 0.f) {
        rap = asinh(pz / xmt);
      } else {
        rap = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      iy = 1 + fem::fint(fem::abs(rap) / dy);
      // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 100
      if (iy < 1 || iy > 50) {
        goto statement_100;
      }
      dyp1(iy) += 1.0f;
      deyp1(iy) += xmt;
      if (ityp == 21) {
        dyg1(iy) += 1.0f;
        deyg1(iy) += xmt;
      }
    statement_100:
      imt = 1 + fem::fint(dxmt / dmt);
      if (rap > ymax || rap <= ymin) {
        goto statement_200;
      }
      if (imt > 200) {
        goto statement_200;
      }
      dmyp1(imt) += 1.0f / xmt;
      if (ityp == 21) {
        dmyg1(imt) += 1.0f / xmt;
      }
    statement_200:;
    }
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      ityp = kftj(i, j);
      px = pjtx(i, j);
      py = pjty(i, j);
      pz = pjtz(i, j);
      pe = pjte(i, j);
      pm = pjtm(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
      dxmt = xmt - pm;
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. PE) THEN
      // C               PRINT *, ' IN HJANA1, TARG STR ', I, ' PART ', J
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', PE
      // C               PRINT *, ' XM = ', PM
      // C               GOTO 400
      // C            END IF
      // C            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      if (xmt > 0.f) {
        rap = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA1 mt=0";
        rap = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      iy = 1 + fem::fint(fem::abs(rap) / dy);
      // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 300
      if (iy < 1 || iy > 50) {
        goto statement_300;
      }
      dyp1(iy) += 1.0f;
      deyp1(iy) += xmt;
      if (ityp == 21) {
        dyg1(iy) += 1.0f;
        deyg1(iy) += xmt;
      }
    statement_300:
      if (rap > ymax || rap <= ymin) {
        goto statement_400;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 200) {
        goto statement_400;
      }
      dmyp1(imt) += 1.0f / xmt;
      if (ityp == 21) {
        dmyg1(imt) += 1.0f / xmt;
      }
    statement_400:;
    }
  }
  // C
  FEM_DO_SAFE(i, 1, nsg) {
    FEM_DO_SAFE(j, 1, njsg(i)) {
      ityp = k2sg(i, j);
      px = pxsg(i, j);
      py = pysg(i, j);
      pz = pzsg(i, j);
      pe = pesg(i, j);
      pm = pmsg(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
      dxmt = xmt - pm;
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. PE) THEN
      // C               PRINT *, ' IN HJANA1, INDP STR ', I, ' PART ', J
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', PE
      // C               PRINT *, ' XM = ', PM
      // C               GOTO 600
      // C            END IF
      // C            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      if (xmt > 0.f) {
        rap = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA1 mt=0";
        rap = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      iy = 1 + fem::fint(fem::abs(rap) / dy);
      // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 500
      if (iy < 1 || iy > 50) {
        goto statement_500;
      }
      dyp1(iy) += 1.0f;
      deyp1(iy) += xmt;
      if (ityp == 21) {
        dyg1(iy) += 1.0f;
        deyg1(iy) += xmt;
      }
    statement_500:
      if (rap > ymax || rap <= ymin) {
        goto statement_600;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 200) {
        goto statement_600;
      }
      dmyp1(imt) += 1.0f / xmt;
      if (ityp == 21) {
        dmyg1(imt) += 1.0f / xmt;
      }
    statement_600:;
    }
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    yr = fem::sqrt(fem::pow2(yp(1, i)) + fem::pow2(yp(2, i)));
    ir = 1 + fem::fint(yr / dr);
    // Clin-4/2008 protect against out-of-bound errors:
    // C         IF (IR .GT. 50) GOTO 601
    if (ir > 50 || ir < 1) {
      goto statement_601;
    }
    dnrpj1(ir) += 1.0f;
    dnrtt1(ir) += 1.0f;
  statement_601:;
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    yr = fem::sqrt(fem::pow2(yt(1, i)) + fem::pow2(yt(2, i)));
    ir = 1 + fem::fint(yr / dr);
    if (ir > 50 || ir < 1) {
      goto statement_602;
    }
    dnrtg1(ir) += 1.0f;
    dnrtt1(ir) += 1.0f;
  statement_602:;
  }
  // C
  FEM_DO_SAFE(i, 1, nsg) {
    y1 = 0.5f * (yp(1, iasg(i, 1)) + yt(1, iasg(i, 2)));
    y2 = 0.5f * (yp(2, iasg(i, 1)) + yt(2, iasg(i, 2)));
    yr = fem::sqrt(fem::pow2(y1) + fem::pow2(y2));
    ir = 1 + fem::fint(yr / dr);
    if (ir > 50 || ir < 1) {
      goto statement_603;
    }
    dnrin1(ir) += 1.0f;
    dnrtt1(ir) += 1.0f;
  statement_603:;
  }
  // C
  FEM_DO_SAFE(i, 1, cmn.mul) {
    ityp = ityp0(i);
    px = fem::sngl(px0(i));
    py = fem::sngl(py0(i));
    pz = fem::sngl(pz0(i));
    pe = fem::sngl(e0(i));
    pm = fem::sngl(xmass0(i));
    xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
    dxmt = xmt - pm;
    // Clin-9/2012 determine rapidity more generally:
    // C         IF (ABS(PZ) .GE. PE) THEN
    // C            PRINT *, ' IN HJANA1, GLUON ', I
    // C            PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
    // C            PRINT *, ' PZ = ', PZ, ' EE = ', PE
    // C            PRINT *, ' XM = ', PM
    // C            GOTO 800
    // C         END IF
    // C         RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
    if (xmt > 0.f) {
      rap = asinh(pz / xmt);
    } else {
      write(6, star), " IN HJANA1 mt=0";
      rap = 1000000.0f * fem::sign(1.f, pz);
    }
    // C
    iy = 1 + fem::fint(fem::abs(rap) / dy);
    // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
    // C         IF (IY .GT. 50) GOTO 700
    if (iy < 1 || iy > 50) {
      goto statement_700;
    }
    dyg1c(iy) += 1.0f;
    deyg1c(iy) += xmt;
  statement_700:
    if (rap > ymax || rap <= ymin) {
      goto statement_800;
    }
    imt = 1 + fem::fint(dxmt / dmt);
    if (imt > 50) {
      goto statement_800;
    }
    dmyg1c(imt) += 1.0f / xmt;
  statement_800:;
  }
  // C.....count number of particles
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      nsubp++;
      if (kfpj(i, j) == 21) {
        nsubg++;
      }
    }
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      nsubp++;
      if (kftj(i, j) == 21) {
        nsubg++;
      }
    }
  }
  // C
  FEM_DO_SAFE(i, 1, nsg) {
    FEM_DO_SAFE(j, 1, njsg(i)) {
      nsubp++;
      if (k2sg(i, j) == 21) {
        nsubg++;
      }
    }
  }
  nisg += nsg;
  if (cmn.iout == 1) {
    // Cbzdbg2/16/99
    // C      PRINT *, ' in HJANA1 '
    // C      PRINT *, ' total number of partons = ', nsubp
    // C      PRINT *, ' total number of gluons = ', nsubg, MUL
    // C      PRINT *, ' number of projectile strings = ', IHNT2(1)
    // C      PRINT *, ' number of target strings = ', IHNT2(3)
    // C      PRINT *, ' number of independent strings = ', NSG
    write(6, star), " in HJANA1 ";
    write(6, star), " total number of partons = ", nsubp / iw;
    write(6, star), " total number of gluons = ", nsubg / iw;
    // C      PRINT *, ' number of projectile strings = ', IHNT2(1)
    // C      PRINT *, ' number of target strings = ', IHNT2(3)
    write(6, star), " number of independent strings = ", nisg / iw;
    // Cbzdbg2/16/99end
  }
  // C
}

struct hjan1b_save {
  float diff2;
  arr<float> dnrg1b;
  arr<float> dtg1b;
  float gx0;
  float gy0;
  int i;
  int ir;
  int isevt;
  int isrun;
  int it;
  int iw;
  int j;
  int k;
  float r0;
  arr<float> snrg1b;
  arr<float> stg1b;
  float tau7;

  hjan1b_save()
      : diff2(fem::float0),
        dnrg1b(dimension(50), fem::fill0),
        dtg1b(dimension(50), fem::fill0),
        gx0(fem::float0),
        gy0(fem::float0),
        i(fem::int0),
        ir(fem::int0),
        isevt(fem::int0),
        isrun(fem::int0),
        it(fem::int0),
        iw(fem::int0),
        j(fem::int0),
        k(fem::int0),
        r0(fem::float0),
        snrg1b(dimension(50), fem::fill0),
        stg1b(dimension(50), fem::fill0),
        tau7(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine in HJAN1A
// C
void hjan1b(common& cmn) {
  FEM_CMN_SVE(hjan1b);
  common_write write(cmn);
  const int maxptn = 400001;
  arr_cref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_cref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_cref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_cref<int> lstrg1(cmn.lstrg1, dimension(maxptn));
  int& nsp = cmn.nsp;
  int& nst = cmn.nst;
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  const int maxstr = 150001;
  arr_cref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  //
  float& diff2 = sve.diff2;
  arr_ref<float> dnrg1b(sve.dnrg1b, dimension(50));
  arr_ref<float> dtg1b(sve.dtg1b, dimension(50));
  float& gx0 = sve.gx0;
  float& gy0 = sve.gy0;
  int& i = sve.i;
  int& ir = sve.ir;
  int& isevt = sve.isevt;
  int& isrun = sve.isrun;
  int& it = sve.it;
  int& iw = sve.iw;
  int& j = sve.j;
  int& k = sve.k;
  float& r0 = sve.r0;
  arr_ref<float> snrg1b(sve.snrg1b, dimension(50));
  arr_ref<float> stg1b(sve.stg1b, dimension(50));
  float& tau7 = sve.tau7;
  if (is_called_first_time) {
    iw = 0;
  }
  const float dr = 0.2f;
  const float dt = 0.2f;
  // Cc      SAVE /PARA1/
  // Cc      SAVE /prec2/
  // Cc      SAVE /ilist8/
  // Cc      SAVE /SREC1/
  // Cc      SAVE /hjcrdn/
  // Cc      SAVE /HJJET2/
  // Cc      SAVE /AREVT/
  // Cc      SAVE /AROUT/
  // C
  if (isevt == iaevt && isrun == iarun) {
    FEM_DO_SAFE(i, 1, 50) {
      dnrg1b(i) = snrg1b(i);
      dtg1b(i) = stg1b(i);
    }
  } else {
    FEM_DO_SAFE(i, 1, 50) {
      snrg1b(i) = dnrg1b(i);
      stg1b(i) = dtg1b(i);
    }
    isevt = iaevt;
    isrun = iarun;
    iw++;
  }
  // C.....analysis
  FEM_DO_SAFE(i, 1, cmn.mul) {
    j = lstrg1(i);
    // C
    if (j <= nsp) {
      k = j;
      gx0 = yp(1, j);
      gy0 = yp(2, j);
    } else if (j <= nsp + nst) {
      k = j - nsp;
      gx0 = yt(1, k);
      gy0 = yt(2, k);
    } else {
      k = j - nsp - nst;
      gx0 = 0.5f * (yp(1, iasg(k, 1)) + yt(1, iasg(k, 2)));
      gy0 = 0.5f * (yp(2, iasg(k, 1)) + yt(2, iasg(k, 2)));
    }
    r0 = fem::sqrt(fem::pow2((fem::sngl(gx5(i)) - gx0)) +
                   fem::pow2((fem::sngl(gy5(i)) - gy0)));
    ir = 1 + fem::fint(r0 / dr);
    if (ir > 50 || ir < 1) {
      goto statement_100;
    }
    dnrg1b(ir) += 1.0f;
  statement_100:
    // Clin-9/2015 to avoid Floating-Point Exception:
    // C         TAU7 = SQRT(sngl(FT5(I) ** 2 - GZ5(I) ** 2))
    diff2 = fem::sngl(fem::pow2(ft5(i)) - fem::pow2(gz5(i)));
    if (diff2 < 0.f) {
      write(6, star), "5:I,ft5,gz5,diff2=", i, ft5(i), gz5(i), diff2;
      tau7 = 1e-6f;
    } else {
      tau7 = fem::sqrt(diff2);
    }
    // C
    it = 1 + fem::fint(tau7 / dt);
    if (it > 50 || it < 1) {
      goto statement_200;
    }
    dtg1b(it) += 1.0f;
  statement_200:;
  }
  // C
}

struct hjan1a_save {
  arr<float> dgxg1a;
  arr<float> dgyg1a;
  float diff2;
  arr<float> dtg1a;
  int i;
  int igx;
  int igy;
  int isevt;
  int isrun;
  int it;
  int iw;
  arr<float> sgxg1a;
  arr<float> sgyg1a;
  arr<float> stg1a;

  hjan1a_save()
      : dgxg1a(dimension(50), fem::fill0),
        dgyg1a(dimension(50), fem::fill0),
        diff2(fem::float0),
        dtg1a(dimension(50), fem::fill0),
        i(fem::int0),
        igx(fem::int0),
        igy(fem::int0),
        isevt(fem::int0),
        isrun(fem::int0),
        it(fem::int0),
        iw(fem::int0),
        sgxg1a(dimension(50), fem::fill0),
        sgyg1a(dimension(50), fem::fill0),
        stg1a(dimension(50), fem::fill0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine in ZPC after generation of additional initial
// C.....phase space distributions.
// C
void hjan1a(common& cmn) {
  FEM_CMN_SVE(hjan1a);
  common_write write(cmn);
  const int maxptn = 400001;
  arr_cref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_cref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_cref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  //
  arr_ref<float> dgxg1a(sve.dgxg1a, dimension(50));
  arr_ref<float> dgyg1a(sve.dgyg1a, dimension(50));
  float& diff2 = sve.diff2;
  arr_ref<float> dtg1a(sve.dtg1a, dimension(50));
  int& i = sve.i;
  int& igx = sve.igx;
  int& igy = sve.igy;
  int& isevt = sve.isevt;
  int& isrun = sve.isrun;
  int& it = sve.it;
  int& iw = sve.iw;
  arr_ref<float> sgxg1a(sve.sgxg1a, dimension(50));
  arr_ref<float> sgyg1a(sve.sgyg1a, dimension(50));
  arr_ref<float> stg1a(sve.stg1a, dimension(50));
  if (is_called_first_time) {
    iw = 0;
  }
  const float dgx = 0.2f;
  const float dgy = 0.2f;
  const float dt = 0.2f;
  // Cc      SAVE /PARA1/
  // Cc      SAVE /prec2/
  // Cc      SAVE /AREVT/
  // Cc      SAVE /AROUT/
  // C
  if (isevt == iaevt && isrun == iarun) {
    FEM_DO_SAFE(i, 1, 50) {
      dgxg1a(i) = sgxg1a(i);
      dgyg1a(i) = sgyg1a(i);
      dtg1a(i) = stg1a(i);
    }
  } else {
    FEM_DO_SAFE(i, 1, 50) {
      sgxg1a(i) = dgxg1a(i);
      sgyg1a(i) = dgyg1a(i);
      stg1a(i) = dtg1a(i);
    }
    isevt = iaevt;
    isrun = iarun;
    iw++;
  }
  // C.....analysis
  FEM_DO_SAFE(i, 1, cmn.mul) {
    igx = 1 + fem::fint(fem::sngl(fem::abs(gx5(i))) / dgx);
    // Clin-4/2008 protect against out-of-bound errors:
    // C         IF (IGX .GT. 50) GOTO 100
    if (igx > 50 || igx < 1) {
      goto statement_100;
    }
    dgxg1a(igx) += 1.0f;
  statement_100:
    igy = 1 + fem::fint(fem::sngl(fem::abs(gy5(i))) / dgy);
    if (igy > 50 || igy < 1) {
      goto statement_200;
    }
    dgyg1a(igy) += 1.0f;
  statement_200:
    // Clin-9/2015 to avoid Floating-Point Exception:
    // C         IT = 1 + int(sngl(SQRT(FT5(I) ** 2 - GZ5(I) ** 2)) / DT)
    diff2 = fem::sngl(fem::pow2(ft5(i)) - fem::pow2(gz5(i)));
    if (diff2 < 0.f) {
      write(6, star), "1:I,ft5,gz5,diff2=", i, ft5(i), gz5(i), diff2;
      it = 1;
    } else {
      it = 1 + fem::fint(fem::sqrt(diff2) / dt);
    }
    // C
    if (it > 50 || it < 1) {
      goto statement_300;
    }
    dtg1a(it) += 1.0f;
  statement_300:;
  }
  hjan1b(cmn);
  // C
}

struct hjan2a_save {
  arr<float> dgxg2a;
  arr<float> dgxp2a;
  arr<float> dgyg2a;
  arr<float> dgyp2a;
  float diff2;
  arr<float> dtg2a;
  arr<float> dtp2a;
  int i;
  int igx;
  int igy;
  int isevt;
  int isrun;
  int it;
  int iw;
  int j;
  arr<float> sgxg2a;
  arr<float> sgxp2a;
  arr<float> sgyg2a;
  arr<float> sgyp2a;
  arr<float> stg2a;
  arr<float> stp2a;

  hjan2a_save()
      : dgxg2a(dimension(50), fem::fill0),
        dgxp2a(dimension(50), fem::fill0),
        dgyg2a(dimension(50), fem::fill0),
        dgyp2a(dimension(50), fem::fill0),
        diff2(fem::float0),
        dtg2a(dimension(50), fem::fill0),
        dtp2a(dimension(50), fem::fill0),
        i(fem::int0),
        igx(fem::int0),
        igy(fem::int0),
        isevt(fem::int0),
        isrun(fem::int0),
        it(fem::int0),
        iw(fem::int0),
        j(fem::int0),
        sgxg2a(dimension(50), fem::fill0),
        sgxp2a(dimension(50), fem::fill0),
        sgyg2a(dimension(50), fem::fill0),
        sgyp2a(dimension(50), fem::fill0),
        stg2a(dimension(50), fem::fill0),
        stp2a(dimension(50), fem::fill0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....subroutine called by HJANA2
void hjan2a(common& cmn) {
  FEM_CMN_SVE(hjan2a);
  common_write write(cmn);
  const int maxptn = 400001;
  arr_cref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_cref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_cref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  arr_cref<int> npj(cmn.npj, dimension(300));
  arr_cref<int, 2> kfpj(cmn.kfpj, dimension(300, 500));
  arr_cref<int> ntj(cmn.ntj, dimension(300));
  arr_cref<int, 2> kftj(cmn.kftj, dimension(300, 500));
  const int maxstr = 150001;
  arr_cref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_cref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  arr_cref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  //
  arr_ref<float> dgxg2a(sve.dgxg2a, dimension(50));
  arr_ref<float> dgxp2a(sve.dgxp2a, dimension(50));
  arr_ref<float> dgyg2a(sve.dgyg2a, dimension(50));
  arr_ref<float> dgyp2a(sve.dgyp2a, dimension(50));
  float& diff2 = sve.diff2;
  arr_ref<float> dtg2a(sve.dtg2a, dimension(50));
  arr_ref<float> dtp2a(sve.dtp2a, dimension(50));
  int& i = sve.i;
  int& igx = sve.igx;
  int& igy = sve.igy;
  int& isevt = sve.isevt;
  int& isrun = sve.isrun;
  int& it = sve.it;
  int& iw = sve.iw;
  int& j = sve.j;
  arr_ref<float> sgxg2a(sve.sgxg2a, dimension(50));
  arr_ref<float> sgxp2a(sve.sgxp2a, dimension(50));
  arr_ref<float> sgyg2a(sve.sgyg2a, dimension(50));
  arr_ref<float> sgyp2a(sve.sgyp2a, dimension(50));
  arr_ref<float> stg2a(sve.stg2a, dimension(50));
  arr_ref<float> stp2a(sve.stp2a, dimension(50));
  if (is_called_first_time) {
    iw = 0;
  }
  const float dgx = 0.2f;
  const float dgy = 0.2f;
  const float dt = 0.2f;
  // C
  // Cc      SAVE /PARA1/
  // Cc      SAVE /prec2/
  // Cc      SAVE /HPARNT/
  // Cc      SAVE /hjcrdn/
  // Cc      SAVE /HJJET1/
  // Cc      SAVE /HJJET2/
  // Cc      SAVE /AREVT/
  // Cc      SAVE /AROUT/
  // C
  if (isevt == iaevt && isrun == iarun) {
    FEM_DO_SAFE(i, 1, 50) {
      dgxp2a(i) = sgxp2a(i);
      dgyp2a(i) = sgyp2a(i);
      dtp2a(i) = stp2a(i);
      dgxg2a(i) = sgxg2a(i);
      dgyg2a(i) = sgyg2a(i);
      dtg2a(i) = stg2a(i);
    }
  } else {
    FEM_DO_SAFE(i, 1, 50) {
      sgxp2a(i) = dgxp2a(i);
      sgyp2a(i) = dgyp2a(i);
      stp2a(i) = dtp2a(i);
      sgxg2a(i) = dgxg2a(i);
      sgyg2a(i) = dgyg2a(i);
      stg2a(i) = dtg2a(i);
    }
    isevt = iaevt;
    isrun = iarun;
    iw++;
  }
  // C.....analysis
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      if (kfpj(i, j) != 21) {
        igx = 1 + fem::fint(fem::abs(yp(1, i)) / dgx);
        if (igx > 50 || igx < 1) {
          goto statement_100;
        }
        dgxp2a(igx) += 1.0f;
      statement_100:
        igy = 1 + fem::fint(fem::abs(yp(2, i)) / dgy);
        if (igy > 50 || igy < 1) {
          goto statement_200;
        }
        dgyp2a(igy) += 1.0f;
      statement_200:
        it = 1;
        dtp2a(it) += 1.0f;
      }
    }
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      if (kftj(i, j) != 21) {
        igx = 1 + fem::fint(fem::abs(yt(1, i)) / dgx);
        if (igx > 50 || igx < 1) {
          goto statement_300;
        }
        dgxp2a(igx) += 1.0f;
      statement_300:
        igy = 1 + fem::fint(fem::abs(yt(2, i)) / dgy);
        if (igy > 50 || igy < 1) {
          goto statement_400;
        }
        dgyp2a(igy) += 1.0f;
      statement_400:
        it = 1;
        dtp2a(it) += 1.0f;
      }
    }
  }
  // C
  FEM_DO_SAFE(i, 1, cmn.nsg) {
    FEM_DO_SAFE(j, 1, njsg(i)) {
      if (k2sg(i, j) != 21) {
        igx = 1 + fem::fint(
                      fem::abs(0.5f * (yp(1, iasg(i, 1)) + yt(1, iasg(i, 2)))) /
                      dgx);
        if (igx > 50 || igx < 1) {
          goto statement_500;
        }
        dgxp2a(igx) += 1.0f;
      statement_500:
        igy = 1 + fem::fint(
                      fem::abs(0.5f * (yp(2, iasg(i, 1)) + yt(2, iasg(i, 2)))) /
                      dgy);
        if (igy > 50 || igy < 1) {
          goto statement_600;
        }
        dgyp2a(igy) += 1.0f;
      statement_600:
        it = 1;
        dtp2a(it) += 1.0f;
      }
    }
  }
  // C
  FEM_DO_SAFE(i, 1, cmn.mul) {
    igx = 1 + fem::fint(fem::abs(fem::sngl(gx5(i))) / dgx);
    if (igx > 50 || igx < 1) {
      goto statement_700;
    }
    dgxg2a(igx) += 1.0f;
    dgxp2a(igx) += 1.0f;
  statement_700:
    igy = 1 + fem::fint(fem::abs(fem::sngl(gy5(i))) / dgy);
    if (igy > 50 || igy < 1) {
      goto statement_800;
    }
    dgyg2a(igy) += 1.0f;
    dgyp2a(igy) += 1.0f;
  statement_800:
    // Clin-9/2015 to avoid Floating-Point Exception:
    // C         IT = 1 + int(SQRT(sngl(FT5(I) ** 2 - GZ5(I) ** 2)) / DT)
    diff2 = fem::sngl(fem::pow2(ft5(i)) - fem::pow2(gz5(i)));
    if (diff2 < 0.f) {
      write(6, star), "3:I,ft5,gz5,diff2=", i, ft5(i), gz5(i), diff2;
      it = 1;
    } else {
      it = 1 + fem::fint(fem::sqrt(diff2) / dt);
    }
    // C
    if (it > 50 || it < 1) {
      goto statement_900;
    }
    dtg2a(it) += 1.0f;
    dtp2a(it) += 1.0f;
  statement_900:;
  }
  // C
}

struct hjan2b_save {
  float diff2;
  arr<float> dnrg2b;
  float dtau;
  arr<float> dtg2b;
  float gx0;
  float gy0;
  int i;
  int ir;
  int isevt;
  int isrun;
  int it;
  int iw;
  int j;
  float r0;
  arr<float> snrg2b;
  arr<float> stg2b;
  float tau7;

  hjan2b_save()
      : diff2(fem::float0),
        dnrg2b(dimension(50), fem::fill0),
        dtau(fem::float0),
        dtg2b(dim1(-24, 25), fem::fill0),
        gx0(fem::float0),
        gy0(fem::float0),
        i(fem::int0),
        ir(fem::int0),
        isevt(fem::int0),
        isrun(fem::int0),
        it(fem::int0),
        iw(fem::int0),
        j(fem::int0),
        r0(fem::float0),
        snrg2b(dimension(50), fem::fill0),
        stg2b(dim1(-24, 25), fem::fill0),
        tau7(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine in HJANA2
// C
void hjan2b(common& cmn) {
  FEM_CMN_SVE(hjan2b);
  common_write write(cmn);
  const int maxptn = 400001;
  arr_cref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_cref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_cref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_cref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_cref<int> lstrg1(cmn.lstrg1, dimension(maxptn));
  const int maxstr = 150001;
  arr_cref<double> ataui(cmn.ataui, dimension(maxstr));
  arr_cref<double> zt1(cmn.zt1, dimension(maxstr));
  arr_cref<double> zt2(cmn.zt2, dimension(maxstr));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  //
  float& diff2 = sve.diff2;
  arr_ref<float> dnrg2b(sve.dnrg2b, dimension(50));
  float& dtau = sve.dtau;
  arr_ref<float> dtg2b(sve.dtg2b, dim1(-24, 25));
  float& gx0 = sve.gx0;
  float& gy0 = sve.gy0;
  int& i = sve.i;
  int& ir = sve.ir;
  int& isevt = sve.isevt;
  int& isrun = sve.isrun;
  int& it = sve.it;
  int& iw = sve.iw;
  int& j = sve.j;
  float& r0 = sve.r0;
  arr_ref<float> snrg2b(sve.snrg2b, dimension(50));
  arr_ref<float> stg2b(sve.stg2b, dim1(-24, 25));
  float& tau7 = sve.tau7;
  if (is_called_first_time) {
    iw = 0;
  }
  const float dr = 0.2f;
  const float dt = 0.2f;
  // C
  // Cc      SAVE /PARA1/
  // Cc      SAVE /prec2/
  // Cc      SAVE /ilist8/
  // Cc      SAVE /SREC1/
  // Cc      SAVE /SREC2/
  // Cc      SAVE /hjcrdn/
  // Cc      SAVE /HJJET2/
  // Cc      SAVE /AREVT/
  // Cc      SAVE /AROUT/
  // C
  if (isevt == iaevt && isrun == iarun) {
    FEM_DO_SAFE(i, 1, 50) {
      dnrg2b(i) = snrg2b(i);
      dtg2b(i - 25) = stg2b(i - 25);
    }
  } else {
    FEM_DO_SAFE(i, 1, 50) {
      snrg2b(i) = dnrg2b(i);
      stg2b(i - 25) = dtg2b(i - 25);
    }
    isevt = iaevt;
    isrun = iarun;
    iw++;
  }
  // C.....analysis
  FEM_DO_SAFE(i, 1, cmn.mul) {
    j = lstrg1(i);
    gx0 = fem::sngl(zt1(j));
    gy0 = fem::sngl(zt2(j));
    r0 = fem::sqrt(fem::pow2((fem::sngl(gx5(i)) - gx0)) +
                   fem::pow2((fem::sngl(gy5(i)) - gy0)));
    ir = 1 + fem::fint(r0 / dr);
    if (ir > 50 || ir < 1) {
      goto statement_100;
    }
    dnrg2b(ir) += 1.0f;
  statement_100:
    // Clin-9/2015 to avoid Floating-Point Exception:
    // C         TAU7 = SQRT(sngl(FT5(I) ** 2 - GZ5(I) ** 2))
    diff2 = fem::sngl(fem::pow2(ft5(i)) - fem::pow2(gz5(i)));
    if (diff2 < 0.f) {
      write(6, star), "4:I,ft5,gz5,diff2=", i, ft5(i), gz5(i), diff2;
      tau7 = 1e-6f;
    } else {
      tau7 = fem::sqrt(diff2);
    }
    // C
    dtau = tau7 - fem::sngl(ataui(j));
    it = 1 + fem::fint(dtau / dt);
    // Cbzdbg2/21/99
    // C         IF (ABS(IT) .GT. 25) GOTO 200
    if (it > 25 || it < -24) {
      goto statement_200;
    }
    // Cbzdbg2/21/99end
    dtg2b(it) += 1.0f;
  statement_200:;
  }
  // C
}

struct hjana2_save {
  arr<float> deyg2;
  arr<float> deyg2c;
  arr<float> deyp2;
  arr<float> dmyg2;
  arr<float> dmyg2c;
  arr<float> dmyp2;
  arr<float> dnrin2;
  arr<float> dnrpj2;
  arr<float> dnrtg2;
  arr<float> dnrtt2;
  arr<float> dtin2;
  arr<float> dtpj2;
  arr<float> dttg2;
  arr<float> dttot2;
  float dxmt;
  arr<float> dyg2;
  arr<float> dyg2c;
  arr<float> dyp2;
  int i;
  int imt;
  int ir;
  int isevt;
  int isrun;
  int it;
  int ityp;
  int iw;
  int iy;
  int j;
  int nisg;
  int nisgs;
  int nj;
  int nsubg;
  int nsubgs;
  int nsubp;
  int nsubps;
  float pe;
  float pm;
  float px;
  float py;
  float pz;
  float rap;
  arr<float> seyg2;
  arr<float> seyg2c;
  arr<float> seyp2;
  arr<float> smyg2;
  arr<float> smyg2c;
  arr<float> smyp2;
  arr<float> snrin2;
  arr<float> snrpj2;
  arr<float> snrtg2;
  arr<float> snrtt2;
  arr<float> snyg2;
  arr<float> snyg2c;
  arr<float> snyp2;
  arr<float> stin2;
  arr<float> stpj2;
  arr<float> sttg2;
  arr<float> sttot2;
  float xmt;
  float yr;

  hjana2_save()
      : deyg2(dimension(50), fem::fill0),
        deyg2c(dimension(50), fem::fill0),
        deyp2(dimension(50), fem::fill0),
        dmyg2(dimension(200), fem::fill0),
        dmyg2c(dimension(50), fem::fill0),
        dmyp2(dimension(200), fem::fill0),
        dnrin2(dimension(50), fem::fill0),
        dnrpj2(dimension(50), fem::fill0),
        dnrtg2(dimension(50), fem::fill0),
        dnrtt2(dimension(50), fem::fill0),
        dtin2(dimension(50), fem::fill0),
        dtpj2(dimension(50), fem::fill0),
        dttg2(dimension(50), fem::fill0),
        dttot2(dimension(50), fem::fill0),
        dxmt(fem::float0),
        dyg2(dimension(50), fem::fill0),
        dyg2c(dimension(50), fem::fill0),
        dyp2(dimension(50), fem::fill0),
        i(fem::int0),
        imt(fem::int0),
        ir(fem::int0),
        isevt(fem::int0),
        isrun(fem::int0),
        it(fem::int0),
        ityp(fem::int0),
        iw(fem::int0),
        iy(fem::int0),
        j(fem::int0),
        nisg(fem::int0),
        nisgs(fem::int0),
        nj(fem::int0),
        nsubg(fem::int0),
        nsubgs(fem::int0),
        nsubp(fem::int0),
        nsubps(fem::int0),
        pe(fem::float0),
        pm(fem::float0),
        px(fem::float0),
        py(fem::float0),
        pz(fem::float0),
        rap(fem::float0),
        seyg2(dimension(50), fem::fill0),
        seyg2c(dimension(50), fem::fill0),
        seyp2(dimension(50), fem::fill0),
        smyg2(dimension(200), fem::fill0),
        smyg2c(dimension(50), fem::fill0),
        smyp2(dimension(200), fem::fill0),
        snrin2(dimension(50), fem::fill0),
        snrpj2(dimension(50), fem::fill0),
        snrtg2(dimension(50), fem::fill0),
        snrtt2(dimension(50), fem::fill0),
        snyg2(dimension(50), fem::fill0),
        snyg2c(dimension(50), fem::fill0),
        snyp2(dimension(50), fem::fill0),
        stin2(dimension(50), fem::fill0),
        stpj2(dimension(50), fem::fill0),
        sttg2(dimension(50), fem::fill0),
        sttot2(dimension(50), fem::fill0),
        xmt(fem::float0),
        yr(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine in HIJING after parton cascade evolution
void hjana2(common& cmn) {
  FEM_CMN_SVE(hjana2);
  common_write write(cmn);
  arr_cref<int> ihnt2(cmn.ihnt2, dimension(50));
  const int maxstr = 150001;
  arr_cref<double> ataui(cmn.ataui, dimension(maxstr));
  arr_cref<double> zt1(cmn.zt1, dimension(maxstr));
  arr_cref<double> zt2(cmn.zt2, dimension(maxstr));
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
  int& nsg = cmn.nsg;
  arr_cref<int> njsg(cmn.njsg, dimension(maxstr));
  arr_cref<int, 2> k2sg(cmn.k2sg, dimension(maxstr, 100));
  arr_cref<float, 2> pxsg(cmn.pxsg, dimension(maxstr, 100));
  arr_cref<float, 2> pysg(cmn.pysg, dimension(maxstr, 100));
  arr_cref<float, 2> pzsg(cmn.pzsg, dimension(maxstr, 100));
  arr_cref<float, 2> pesg(cmn.pesg, dimension(maxstr, 100));
  arr_cref<float, 2> pmsg(cmn.pmsg, dimension(maxstr, 100));
  const int maxptn = 400001;
  arr_cref<double> px5(cmn.px5, dimension(maxptn));
  arr_cref<double> py5(cmn.py5, dimension(maxptn));
  arr_cref<double> pz5(cmn.pz5, dimension(maxptn));
  arr_cref<double> e5(cmn.e5, dimension(maxptn));
  arr_cref<double> xmass5(cmn.xmass5, dimension(maxptn));
  arr_cref<int> ityp5(cmn.ityp5, dimension(maxptn));
  int& iaevt = cmn.iaevt;
  int& iarun = cmn.iarun;
  int& isoft = cmn.isoft;
  arr_cref<double, 2> pxsgs(cmn.pxsgs, dimension(maxstr, 3));
  arr_cref<double, 2> pysgs(cmn.pysgs, dimension(maxstr, 3));
  arr_cref<double, 2> pzsgs(cmn.pzsgs, dimension(maxstr, 3));
  arr_cref<double, 2> pesgs(cmn.pesgs, dimension(maxstr, 3));
  arr_cref<double, 2> pmsgs(cmn.pmsgs, dimension(maxstr, 3));
  arr_cref<int, 2> k2sgs(cmn.k2sgs, dimension(maxstr, 3));
  arr_cref<int> njsgs(cmn.njsgs, dimension(maxstr));
  //
  arr_ref<float> deyg2(sve.deyg2, dimension(50));
  arr_ref<float> deyg2c(sve.deyg2c, dimension(50));
  arr_ref<float> deyp2(sve.deyp2, dimension(50));
  arr_ref<float> dmyg2(sve.dmyg2, dimension(200));
  arr_ref<float> dmyg2c(sve.dmyg2c, dimension(50));
  arr_ref<float> dmyp2(sve.dmyp2, dimension(200));
  arr_ref<float> dnrin2(sve.dnrin2, dimension(50));
  arr_ref<float> dnrpj2(sve.dnrpj2, dimension(50));
  arr_ref<float> dnrtg2(sve.dnrtg2, dimension(50));
  arr_ref<float> dnrtt2(sve.dnrtt2, dimension(50));
  arr_ref<float> dtin2(sve.dtin2, dimension(50));
  arr_ref<float> dtpj2(sve.dtpj2, dimension(50));
  arr_ref<float> dttg2(sve.dttg2, dimension(50));
  arr_ref<float> dttot2(sve.dttot2, dimension(50));
  float& dxmt = sve.dxmt;
  arr_ref<float> dyg2(sve.dyg2, dimension(50));
  arr_ref<float> dyg2c(sve.dyg2c, dimension(50));
  arr_ref<float> dyp2(sve.dyp2, dimension(50));
  int& i = sve.i;
  int& imt = sve.imt;
  int& ir = sve.ir;
  int& isevt = sve.isevt;
  int& isrun = sve.isrun;
  int& it = sve.it;
  int& ityp = sve.ityp;
  int& iw = sve.iw;
  int& iy = sve.iy;
  int& j = sve.j;
  int& nisg = sve.nisg;
  int& nisgs = sve.nisgs;
  int& nj = sve.nj;
  int& nsubg = sve.nsubg;
  int& nsubgs = sve.nsubgs;
  int& nsubp = sve.nsubp;
  int& nsubps = sve.nsubps;
  float& pe = sve.pe;
  float& pm = sve.pm;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& rap = sve.rap;
  arr_ref<float> seyg2(sve.seyg2, dimension(50));
  arr_ref<float> seyg2c(sve.seyg2c, dimension(50));
  arr_ref<float> seyp2(sve.seyp2, dimension(50));
  arr_ref<float> smyg2(sve.smyg2, dimension(200));
  arr_ref<float> smyg2c(sve.smyg2c, dimension(50));
  arr_ref<float> smyp2(sve.smyp2, dimension(200));
  arr_ref<float> snrin2(sve.snrin2, dimension(50));
  arr_ref<float> snrpj2(sve.snrpj2, dimension(50));
  arr_ref<float> snrtg2(sve.snrtg2, dimension(50));
  arr_ref<float> snrtt2(sve.snrtt2, dimension(50));
  arr_ref<float> snyg2(sve.snyg2, dimension(50));
  arr_ref<float> snyg2c(sve.snyg2c, dimension(50));
  arr_ref<float> snyp2(sve.snyp2, dimension(50));
  arr_ref<float> stin2(sve.stin2, dimension(50));
  arr_ref<float> stpj2(sve.stpj2, dimension(50));
  arr_ref<float> sttg2(sve.sttg2, dimension(50));
  arr_ref<float> sttot2(sve.sttot2, dimension(50));
  float& xmt = sve.xmt;
  float& yr = sve.yr;
  if (is_called_first_time) {
    iw = 0;
  }
  const float dy = 0.2f;
  const float ymax = 1.0f;
  const float ymin = -1.0f;
  const float dmt = 0.05f;
  const float dr = 0.2f;
  const float dt = 0.2f;
  // C
  // Cc      SAVE /PARA1/
  // Cc      SAVE /HPARNT/
  // Cc      SAVE /SREC2/
  // Cc      SAVE /HJJET1/
  // Cc      SAVE /HJJET2/
  // Cc      SAVE /prec2/
  // Cc      SAVE /AREVT/
  // Cc      SAVE /AROUT/
  // Cc      SAVE /anim/
  // Cc      SAVE /SOFT/
  // C
  if (isevt == iaevt && isrun == iarun) {
    FEM_DO_SAFE(i, 1, 200) {
      dmyp2(i) = smyp2(i);
      dmyg2(i) = smyg2(i);
    }
    // C
    FEM_DO_SAFE(i, 1, 50) {
      dyp2(i) = snyp2(i);
      deyp2(i) = seyp2(i);
      dyg2(i) = snyg2(i);
      deyg2(i) = seyg2(i);
      dnrpj2(i) = snrpj2(i);
      dnrtg2(i) = snrtg2(i);
      dnrin2(i) = snrin2(i);
      dnrtt2(i) = snrtt2(i);
      dtpj2(i) = stpj2(i);
      dttg2(i) = sttg2(i);
      dtin2(i) = stin2(i);
      dttot2(i) = sttot2(i);
      dyg2c(i) = snyg2c(i);
      dmyg2c(i) = smyg2c(i);
      deyg2c(i) = seyg2c(i);
    }
    nsubp = nsubps;
    nsubg = nsubgs;
    nisg = nisgs;
  } else {
    FEM_DO_SAFE(i, 1, 200) {
      smyp2(i) = dmyp2(i);
      smyg2(i) = dmyg2(i);
    }
    // C
    FEM_DO_SAFE(i, 1, 50) {
      snyp2(i) = dyp2(i);
      seyp2(i) = deyp2(i);
      snyg2(i) = dyg2(i);
      seyg2(i) = deyg2(i);
      snrpj2(i) = dnrpj2(i);
      snrtg2(i) = dnrtg2(i);
      snrin2(i) = dnrin2(i);
      snrtt2(i) = dnrtt2(i);
      stpj2(i) = dtpj2(i);
      sttg2(i) = dttg2(i);
      stin2(i) = dtin2(i);
      sttot2(i) = dttot2(i);
      snyg2c(i) = dyg2c(i);
      smyg2c(i) = dmyg2c(i);
      seyg2c(i) = deyg2c(i);
    }
    nsubps = nsubp;
    nsubgs = nsubg;
    nisgs = nisg;
    isevt = iaevt;
    isrun = iarun;
    iw++;
  }
  // C
  // Clin-4/28/01:
  if (isoft == 3 || isoft == 4 || isoft == 5) {
    goto statement_510;
  }
  // C
  // C.....analysis
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      ityp = kfpj(i, j);
      px = pjpx(i, j);
      py = pjpy(i, j);
      pz = pjpz(i, j);
      pe = pjpe(i, j);
      pm = pjpm(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
      dxmt = xmt - pm;
      // Clin-9/2012 determine rapidity more generally:
      // Ccbzdbg2/16/99
      // Cc            IF (ABS(PZ) .GE. PE) GOTO 200
      // C            IF (ABS(PZ) .GE. PE) THEN
      // C               PRINT *, ' IN HJANA2, PROJ STR ', I, ' PART ', J
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', PE
      // C               PRINT *, ' XM = ', PM
      // C               GOTO 200
      // C            END IF
      // Ccbzdbg2/16/99end
      // C            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      if (xmt > 0.f) {
        rap = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA2 mt=0";
        rap = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      iy = 1 + fem::fint(fem::abs(rap) / dy);
      // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 100
      if (iy < 1 || iy > 50) {
        goto statement_100;
      }
      dyp2(iy) += 1.0f;
      deyp2(iy) += xmt;
      if (ityp == 21) {
        dyg2(iy) += 1.0f;
        deyg2(iy) += xmt;
      }
    statement_100:
      if (rap > ymax || rap <= ymin) {
        goto statement_200;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 200) {
        goto statement_200;
      }
      dmyp2(imt) += 1.0f / xmt;
      if (ityp == 21) {
        dmyg2(imt) += 1.0f / xmt;
      }
    statement_200:;
    }
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      ityp = kftj(i, j);
      px = pjtx(i, j);
      py = pjty(i, j);
      pz = pjtz(i, j);
      pe = pjte(i, j);
      pm = pjtm(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
      dxmt = xmt - pm;
      // Clin-9/2012 determine rapidity more generally:
      // Ccbzdbg2/16/99
      // Cc            IF (ABS(PZ) .GE. PE) GOTO 400
      // C            IF (ABS(PZ) .GE. PE) THEN
      // C               PRINT *, ' IN HJANA2, TARG STR ', I, ' PART ', J
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', PE
      // C               PRINT *, ' XM = ', PM
      // C               GOTO 400
      // C            END IF
      // Ccbzdbg2/16/99end
      // C            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      if (xmt > 0.f) {
        rap = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA2 mt=0";
        rap = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      iy = 1 + fem::fint(fem::abs(rap) / dy);
      // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 300
      if (iy < 1 || iy > 50) {
        goto statement_300;
      }
      dyp2(iy) += 1.0f;
      deyp2(iy) += xmt;
      if (ityp == 21) {
        dyg2(iy) += 1.0f;
        deyg2(iy) += xmt;
      }
    statement_300:
      if (rap > ymax || rap <= ymin) {
        goto statement_400;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 200) {
        goto statement_400;
      }
      dmyp2(imt) += 1.0f / xmt;
      if (ityp == 21) {
        dmyg2(imt) += 1.0f / xmt;
      }
    statement_400:;
    }
  }
// C
// Clin-4/28/01:
statement_510:
  // C
  FEM_DO_SAFE(i, 1, nsg) {
    // Clin-4/25/01 soft3:
    // C         DO J = 1, NJSG(I)
    nj = njsg(i);
    if (isoft == 3 || isoft == 4 || isoft == 5) {
      nj = njsgs(i);
    }
    FEM_DO_SAFE(j, 1, nj) {
      // Clin-4/25/01-end
      // C
      ityp = k2sg(i, j);
      px = pxsg(i, j);
      py = pysg(i, j);
      pz = pzsg(i, j);
      pe = pesg(i, j);
      pm = pmsg(i, j);
      // Clin-4/25/01 soft3:
      if (isoft == 3 || isoft == 4 || isoft == 5) {
        ityp = k2sgs(i, j);
        px = fem::sngl(pxsgs(i, j));
        py = fem::sngl(pysgs(i, j));
        pz = fem::sngl(pzsgs(i, j));
        pe = fem::sngl(pesgs(i, j));
        pm = fem::sngl(pmsgs(i, j));
      }
      // Clin-4/25/01-end
      // C
      // Clin-9/2012 determine rapidity more generally:
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
      dxmt = xmt - pm;
      // Ccbzdbg2/16/99
      // Cc            IF (ABS(PZ) .GE. PE) GOTO 600
      // C            IF (ABS(PZ) .GE. PE) THEN
      // C               PRINT *, ' IN HJANA2, INDP STR ', I, ' PART ', J
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', PE
      // C               PRINT *, ' XM = ', PM
      // C               GOTO 600
      // C            END IF
      // Ccbzdbg2/16/99end
      // C            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      if (xmt > 0.f) {
        rap = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA2 mt=0";
        rap = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      iy = 1 + fem::fint(fem::abs(rap) / dy);
      // Clin-8/2014 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 500
      if (iy < 1 || iy > 50) {
        goto statement_500;
      }
      dyp2(iy) += 1.0f;
      deyp2(iy) += xmt;
      if (ityp == 21) {
        dyg2(iy) += 1.0f;
        deyg2(iy) += xmt;
      }
    statement_500:
      if (rap > ymax || rap <= ymin) {
        goto statement_600;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 200) {
        goto statement_600;
      }
      dmyp2(imt) += 1.0f / xmt;
      if (ityp == 21) {
        dmyg2(imt) += 1.0f / xmt;
      }
    statement_600:;
    }
  }
  // C
  // Clin-4/28/01:
  if (isoft == 3 || isoft == 4 || isoft == 5) {
    goto statement_520;
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    j = i;
    yr = fem::sqrt(fem::sngl(fem::pow2(zt1(j)) + fem::pow2(zt2(j))));
    ir = 1 + fem::fint(yr / dr);
    if (ir > 50 || ir < 1) {
      goto statement_601;
    }
    dnrpj2(ir) += 1.0f;
    dnrtt2(ir) += 1.0f;
  statement_601:
    it = 1 + fem::fint(fem::sngl(ataui(j)) / dt);
    if (it > 50 || it < 1) {
      goto statement_602;
    }
    dtpj2(it) += 1.0f;
    dttot2(it) += 1.0f;
  statement_602:;
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    j = i + ihnt2(1);
    yr = fem::sqrt(fem::sngl(fem::pow2(zt1(j)) + fem::pow2(zt2(j))));
    ir = 1 + fem::fint(yr / dr);
    if (ir > 50 || ir < 1) {
      goto statement_603;
    }
    dnrtg2(ir) += 1.0f;
    dnrtt2(ir) += 1.0f;
  statement_603:
    it = 1 + fem::fint(fem::sngl(ataui(j)) / dt);
    if (it > 50 || it < 1) {
      goto statement_604;
    }
    dttg2(it) += 1.0f;
    dttot2(it) += 1.0f;
  statement_604:;
  }
// C
// Clin-4/28/01:
statement_520:
  // C
  FEM_DO_SAFE(i, 1, nsg) {
    j = i + ihnt2(1) + ihnt2(3);
    // Clin-4/28/01:
    if (isoft == 3 || isoft == 4 || isoft == 5) {
      j = i;
    }
    // C
    yr = fem::sqrt(fem::sngl(fem::pow2(zt1(j)) + fem::pow2(zt2(j))));
    ir = 1 + fem::fint(yr / dr);
    if (ir > 50 || ir < 1) {
      goto statement_605;
    }
    dnrin2(ir) += 1.0f;
    dnrtt2(ir) += 1.0f;
  statement_605:
    it = 1 + fem::fint(fem::sngl(ataui(j)) / dt);
    if (it > 50 || it < 1) {
      goto statement_606;
    }
    dtin2(it) += 1.0f;
    dttot2(it) += 1.0f;
  statement_606:;
  }
  // C
  FEM_DO_SAFE(i, 1, cmn.mul) {
    ityp = ityp5(i);
    px = fem::sngl(px5(i));
    py = fem::sngl(py5(i));
    pz = fem::sngl(pz5(i));
    pe = fem::sngl(e5(i));
    pm = fem::sngl(xmass5(i));
    // Clin-9/2012 determine rapidity more generally:
    xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(pm));
    dxmt = xmt - pm;
    // Ccbzdbg2/16/99
    // Cc            IF (ABS(PZ) .GE. PE) GOTO 800
    // C         IF (ABS(PZ) .GE. PE) THEN
    // C            PRINT *, ' IN HJANA2, GLUON ', I
    // C            PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
    // C            PRINT *, ' PZ = ', PZ, ' EE = ', PE
    // C            PRINT *, ' XM = ', PM
    // C            GOTO 800
    // C         END IF
    // Ccbzdbg2/16/99end
    // C         RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
    if (xmt > 0.f) {
      rap = asinh(pz / xmt);
    } else {
      write(6, star), " IN HJANA2 mt=0";
      rap = 1000000.0f * fem::sign(1.f, pz);
    }
    // C
    iy = 1 + fem::fint(fem::abs(rap) / dy);
    // Clin-9/2012 prevent possible segmentation fault (due to IY<=0):
    // C         IF (IY .GT. 50) GOTO 700
    if (iy < 1 || iy > 50) {
      goto statement_700;
    }
    dyg2c(iy) += 1.0f;
    deyg2c(iy) += xmt;
  statement_700:
    if (rap > ymax || rap <= ymin) {
      goto statement_800;
    }
    imt = 1 + fem::fint(dxmt / dmt);
    if (imt > 50) {
      goto statement_800;
    }
    dmyg2c(imt) += 1.0f / xmt;
  statement_800:;
  }
  // C
  // Clin-4/25/01 soft3:
  if (isoft == 3 || isoft == 4 || isoft == 5) {
    goto statement_530;
  }
  // C
  // C.....count number of particles
  FEM_DO_SAFE(i, 1, ihnt2(1)) {
    FEM_DO_SAFE(j, 1, npj(i)) {
      nsubp++;
      if (kfpj(i, j) == 21) {
        nsubg++;
      }
    }
  }
  // C
  FEM_DO_SAFE(i, 1, ihnt2(3)) {
    FEM_DO_SAFE(j, 1, ntj(i)) {
      nsubp++;
      if (kftj(i, j) == 21) {
        nsubg++;
      }
    }
  }
// C
// Clin-4/25/01 soft3:
statement_530:
  // C
  FEM_DO_SAFE(i, 1, nsg) {
    // Clin-4/25/01 soft3:
    // C         DO J = 1, NJSG(I)
    nj = njsg(i);
    if (isoft == 3 || isoft == 4 || isoft == 5) {
      nj = njsgs(i);
    }
    FEM_DO_SAFE(j, 1, nj) {
      // Clin-4/25/01-end
      // C
      nsubp++;
      // C
      // Clin-4/25/01
      // C            IF (K2SG(I, J) .EQ. 21) nsubg = nsubg + 1
      if (isoft == 3 || isoft == 4 || isoft == 5) {
        if (k2sgs(i, j) == 21) {
          nsubg++;
        }
      } else {
        if (k2sg(i, j) == 21) {
          nsubg++;
        }
      }
      // Clin-4/25/01-end
    }
  }
  // Cbzdbg2/16/99
  nisg += nsg;
  // C
  if (cmn.iout == 1) {
    // Cbzdbg2/16/99end
    // Cbzdbg2/16/99
    // C      PRINT *, ' in HJANA2 '
    // C      PRINT *, ' total number of partons = ', nsubp
    // C      PRINT *, ' total number of gluons = ', nsubg, MUL
    // C      PRINT *, ' number of projectile strings = ', IHNT2(1)
    // C      PRINT *, ' number of target strings = ', IHNT2(3)
    // C      PRINT *, ' number of independent strings = ', NSG
    write(6, star), " in HJANA2 ";
    write(6, star), " total number of partons = ", nsubp / iw;
    write(6, star), " total number of gluons = ", nsubg / iw;
    // C      PRINT *, ' number of projectile strings = ', IHNT2(1)
    // C      PRINT *, ' number of target strings = ', IHNT2(3)
    write(6, star), " number of independent strings = ", nisg / iw;
  }
  // C
  hjan2a(cmn);
  hjan2b(cmn);
  // C
}

struct hjana3_save {
  arr<float> deyh3;
  arr<float> dmyh3;
  arr<float> dndyh3;
  float dxmt;
  float ee;
  int i;
  int imt;
  int ityp;
  int iw;
  int iy;
  int j;
  float px;
  float py;
  float pz;
  float xm;
  float xmt;
  float y;

  hjana3_save()
      : deyh3(dimension(50), fem::fill0),
        dmyh3(dimension(50), fem::fill0),
        dndyh3(dimension(50), fem::fill0),
        dxmt(fem::float0),
        ee(fem::float0),
        i(fem::int0),
        imt(fem::int0),
        ityp(fem::int0),
        iw(fem::int0),
        iy(fem::int0),
        j(fem::int0),
        px(fem::float0),
        py(fem::float0),
        pz(fem::float0),
        xm(fem::float0),
        xmt(fem::float0),
        y(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine before ARTMN
void hjana3(common& cmn) {
  FEM_CMN_SVE(hjana3);
  common_write write(cmn);
  const int maxr = 1;
  arr_cref<int> multi1(cmn.multi1, dimension(maxr));
  const int maxstr = 150001;
  arr_cref<int, 2> ityp1(cmn.ityp1, dimension(maxstr, maxr));
  arr_cref<float, 2> px1(cmn.px1, dimension(maxstr, maxr));
  arr_cref<float, 2> py1(cmn.py1, dimension(maxstr, maxr));
  arr_cref<float, 2> pz1(cmn.pz1, dimension(maxstr, maxr));
  arr_cref<float, 2> ee1(cmn.ee1, dimension(maxstr, maxr));
  arr_cref<float, 2> xm1(cmn.xm1, dimension(maxstr, maxr));
  //
  arr_ref<float> deyh3(sve.deyh3, dimension(50));
  arr_ref<float> dmyh3(sve.dmyh3, dimension(50));
  arr_ref<float> dndyh3(sve.dndyh3, dimension(50));
  float& dxmt = sve.dxmt;
  int& i = sve.i;
  int& imt = sve.imt;
  int& ityp = sve.ityp;
  int& iw = sve.iw;
  int& iy = sve.iy;
  int& j = sve.j;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& xm = sve.xm;
  float& xmt = sve.xmt;
  float& y = sve.y;
  if (is_called_first_time) {
    iw = 0;
  }
  const float dy = 0.2f;
  const float ymin = -1.0f;
  const float ymax = 1.0f;
  const float dmt = 0.05f;
  // C
  // C.....y cut for mt spectrum
  // Cbz11/7/99 end
  // C.....bin width for mt spectrum and y spectrum
  // Cc      SAVE /RUN/
  // Cc      SAVE /ARERC1/
  // Cc      SAVE /ARPRC1/
  // Cc      SAVE /AROUT/
  // Cc      SAVE /iflow/
  // C
  iw++;
  FEM_DO_SAFE(j, 1, cmn.num) {
    FEM_DO_SAFE(i, 1, multi1(j)) {
      ityp = ityp1(i, j);
      if (ityp > -100 && ityp < 100) {
        goto statement_200;
      }
      px = px1(i, j);
      py = py1(i, j);
      pz = pz1(i, j);
      sve.ee = ee1(i, j);
      xm = xm1(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(xm));
      dxmt = xmt - xm;
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. EE) THEN
      // C               PRINT *, 'IN HJANA3'
      // C               PRINT *, ' PARTICLE ', I, ' RUN ', J, 'PREC ERR'
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', EE
      // C               PRINT *, ' XM = ', XM
      // C               GOTO 200
      // C            END IF
      // C            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
      if (xmt > 0.f) {
        y = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA3 mt=0";
        y = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      // C.....rapidity cut for the rapidity distribution
      // C            IY = 1 + int(ABS(Y) / DY)
      // Clin-8/2014 no rapidity shift here:
      // C            IY = 1 + int((Y+10.) / DY)
      iy = 1 + fem::fint(y / dy);
      // Clin-9/2012 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 100
      if (iy < 1 || iy > 50) {
        goto statement_100;
      }
      dndyh3(iy) += 1.0f;
      deyh3(iy) += xmt;
    statement_100:
      // C.....insert rapidity cut for mt spectrum here
      if (y < ymin || y >= ymax) {
        goto statement_200;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 50) {
        goto statement_200;
      }
      dmyh3(imt) += 1.0f / xmt;
    statement_200:;
    }
  }
  // C
}

struct hjana4_save {
  arr<float> deyh4;
  arr<float> dmyh4;
  arr<float> dndyh4;
  float dxmt;
  float ee;
  int i;
  int imt;
  int ityp;
  int iw;
  int iy;
  int j;
  float px;
  float py;
  float pz;
  float xm;
  float xmt;
  float y;

  hjana4_save()
      : deyh4(dimension(50), fem::fill0),
        dmyh4(dimension(50), fem::fill0),
        dndyh4(dimension(50), fem::fill0),
        dxmt(fem::float0),
        ee(fem::float0),
        i(fem::int0),
        imt(fem::int0),
        ityp(fem::int0),
        iw(fem::int0),
        iy(fem::int0),
        j(fem::int0),
        px(fem::float0),
        py(fem::float0),
        pz(fem::float0),
        xm(fem::float0),
        xmt(fem::float0),
        y(fem::float0) {}
};

// C
// C-----------------------------------------------------------------------
// C
// C.....analysis subroutine after ARTMN
void hjana4(common& cmn) {
  FEM_CMN_SVE(hjana4);
  common_write write(cmn);
  const int maxr = 1;
  arr_cref<int> multi1(cmn.multi1, dimension(maxr));
  const int maxstr = 150001;
  arr_cref<int, 2> ityp1(cmn.ityp1, dimension(maxstr, maxr));
  arr_cref<float, 2> px1(cmn.px1, dimension(maxstr, maxr));
  arr_cref<float, 2> py1(cmn.py1, dimension(maxstr, maxr));
  arr_cref<float, 2> pz1(cmn.pz1, dimension(maxstr, maxr));
  arr_cref<float, 2> ee1(cmn.ee1, dimension(maxstr, maxr));
  arr_cref<float, 2> xm1(cmn.xm1, dimension(maxstr, maxr));
  //
  arr_ref<float> deyh4(sve.deyh4, dimension(50));
  arr_ref<float> dmyh4(sve.dmyh4, dimension(50));
  arr_ref<float> dndyh4(sve.dndyh4, dimension(50));
  float& dxmt = sve.dxmt;
  int& i = sve.i;
  int& imt = sve.imt;
  int& ityp = sve.ityp;
  int& iw = sve.iw;
  int& iy = sve.iy;
  int& j = sve.j;
  float& px = sve.px;
  float& py = sve.py;
  float& pz = sve.pz;
  float& xm = sve.xm;
  float& xmt = sve.xmt;
  float& y = sve.y;
  if (is_called_first_time) {
    iw = 0;
  }
  const float dy = 0.2f;
  const float ymin = -1.0f;
  const float ymax = 1.0f;
  const float dmt = 0.05f;
  // C.....y cut for mt spectrum
  // Cbz11/7/99
  // C      PARAMETER (YMIN = -0.5, YMAX = 0.5)
  // Cbz11/7/99 end
  // C.....bin width for mt spectrum and y spectrum
  // C
  // Cc      SAVE /RUN/
  // Cc      SAVE /ARERC1/
  // Cc      SAVE /ARPRC1/
  // Cc      SAVE /AROUT/
  // Cc      SAVE /fflow/
  // C
  iw++;
  FEM_DO_SAFE(j, 1, cmn.num) {
    FEM_DO_SAFE(i, 1, multi1(j)) {
      ityp = ityp1(i, j);
      if (ityp > -100 && ityp < 100) {
        goto statement_200;
      }
      px = px1(i, j);
      py = py1(i, j);
      pz = pz1(i, j);
      sve.ee = ee1(i, j);
      xm = xm1(i, j);
      xmt = fem::sqrt(fem::pow2(px) + fem::pow2(py) + fem::pow2(xm));
      dxmt = xmt - xm;
      // Clin-9/2012 determine rapidity more generally:
      // C            IF (ABS(PZ) .GE. EE) THEN
      // C               PRINT *, 'IN HJANA4'
      // C               PRINT *, ' PARTICLE ', I, ' RUN ', J, 'PREC ERR'
      // C               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
      // C               PRINT *, ' PZ = ', PZ, ' EE = ', EE
      // C               PRINT *, ' XM = ', XM
      // C               GOTO 200
      // C            END IF
      // C            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
      if (xmt > 0.f) {
        y = asinh(pz / xmt);
      } else {
        write(6, star), " IN HJANA4 mt=0";
        y = 1000000.0f * fem::sign(1.f, pz);
      }
      // C
      // C.....rapidity cut for the rapidity distribution
      // C            IY = 1 + int(ABS(Y) / DY)
      // Clin-8/2014 no rapidity shift here:
      // C            IY = 1 + int((Y+10.) / DY)
      iy = 1 + fem::fint(y / dy);
      // Clin-9/2012 prevent possible segmentation fault (due to IY<=0):
      // C            IF (IY .GT. 50) GOTO 100
      if (iy < 1 || iy > 50) {
        goto statement_100;
      }
      dndyh4(iy) += 1.0f;
      deyh4(iy) += xmt;
    statement_100:
      // C.....insert rapidity cut for mt spectrum here
      if (y < ymin || y >= ymax) {
        goto statement_200;
      }
      imt = 1 + fem::fint(dxmt / dmt);
      if (imt > 50) {
        goto statement_200;
      }
      dmyh4(imt) += 1.0f / xmt;
    statement_200:;
    }
  }
  // C
}

struct zpstrg_save {
  float bb;
  double diff2;
  int i;
  int istrg;
  int j;
  int nstr;
  double shift;
  double tau7;

  zpstrg_save()
      : bb(fem::float0),
        diff2(fem::double0),
        i(fem::int0),
        istrg(fem::int0),
        j(fem::int0),
        nstr(fem::int0),
        shift(fem::double0),
        tau7(fem::double0) {}
};

// C
// C=======================================================================
// C
// C.....subroutine to get average values for different strings
// C
void zpstrg(common& cmn) {
  FEM_CMN_SVE(zpstrg);
  common_write write(cmn);
  // COMMON para1
  int& mul = cmn.mul;
  // COMMON prec2
  const int maxptn = 400001;
  arr_ref<double> gx5(cmn.gx5, dimension(maxptn));
  arr_ref<double> gy5(cmn.gy5, dimension(maxptn));
  arr_ref<double> gz5(cmn.gz5, dimension(maxptn));
  arr_ref<double> ft5(cmn.ft5, dimension(maxptn));
  arr_ref<double> px5(cmn.px5, dimension(maxptn));
  arr_ref<double> py5(cmn.py5, dimension(maxptn));
  arr_ref<double> pz5(cmn.pz5, dimension(maxptn));
  arr_ref<double> e5(cmn.e5, dimension(maxptn));
  arr_ref<double> xmass5(cmn.xmass5, dimension(maxptn));
  arr_ref<int> ityp5(cmn.ityp5, dimension(maxptn));
  // COMMON ilist8
  arr_cref<int> lstrg1(cmn.lstrg1, dimension(maxptn));
  // COMMON srec1
  int& nsp = cmn.nsp;
  int& nst = cmn.nst;
  // COMMON srec2
  const int maxstr = 150001;
  arr_ref<double> ataui(cmn.ataui, dimension(maxstr));
  arr_ref<double> zt1(cmn.zt1, dimension(maxstr));
  arr_ref<double> zt2(cmn.zt2, dimension(maxstr));
  arr_ref<double> zt3(cmn.zt3, dimension(maxstr));
  // COMMON hjcrdn
  arr_cref<float, 2> yp(cmn.yp, dimension(3, 300));
  arr_cref<float, 2> yt(cmn.yt, dimension(3, 300));
  // COMMON hjjet2
  arr_cref<int, 2> iasg(cmn.iasg, dimension(maxstr, 3));
  // COMMON hparnt
  arr_cref<float> hint1(cmn.hint1, dimension(100));
  // COMMON anim
  int& isoft = cmn.isoft;
  // COMMON strg
  arr_ref<int> np(cmn.np, dimension(maxstr));
  // COMMON frzprc
  arr_cref<double> gxfrz(cmn.gxfrz, dimension(maxptn));
  arr_cref<double> gyfrz(cmn.gyfrz, dimension(maxptn));
  arr_cref<double> gzfrz(cmn.gzfrz, dimension(maxptn));
  arr_cref<double> ftfrz(cmn.ftfrz, dimension(maxptn));
  arr_cref<double> pxfrz(cmn.pxfrz, dimension(maxptn));
  arr_cref<double> pyfrz(cmn.pyfrz, dimension(maxptn));
  arr_cref<double> pzfrz(cmn.pzfrz, dimension(maxptn));
  arr_cref<double> efrz(cmn.efrz, dimension(maxptn));
  arr_cref<double> xmfrz(cmn.xmfrz, dimension(maxptn));
  arr_cref<int> idfrz(cmn.idfrz, dimension(maxptn));
  //
  // SAVE
  float& bb = sve.bb;
  double& diff2 = sve.diff2;
  int& i = sve.i;
  int& istrg = sve.istrg;
  int& j = sve.j;
  int& nstr = sve.nstr;
  double& shift = sve.shift;
  double& tau7 = sve.tau7;
  //
  // C
  // C      REAL*4 YP, YT, PXSG, PYSG, PZSG, PESG, PMSG, HIPR1, HINT1, BB
  // C
  // Cc      SAVE /PARA1/
  // Cc      SAVE /prec2/
  // Cc      SAVE /ilist8/
  // Cc      SAVE /SREC1/
  // Cc      SAVE /SREC2/
  // Cc      SAVE /hjcrdn/
  // Cc      SAVE /HJJET2/
  // Cbz6/28/99 flow1
  // Cc      SAVE /HPARNT/
  // Cbz6/28/99 flow1 end
  // Cc      SAVE /anim/
  // Cc      SAVE /strg/
  // Clin-6/06/02 test local freezeout:
  // Cc      SAVE /frzprc/
  // C
  // Clin-6/06/02 test local freezeout for string melting,
  // C     use space-time values at local freezeout saved in /frzprc/:
  if (isoft == 5) {
    FEM_DO_SAFE(i, 1, mul) {
      ityp5(i) = idfrz(i);
      gx5(i) = gxfrz(i);
      gy5(i) = gyfrz(i);
      gz5(i) = gzfrz(i);
      ft5(i) = ftfrz(i);
      px5(i) = pxfrz(i);
      py5(i) = pyfrz(i);
      pz5(i) = pzfrz(i);
      e5(i) = efrz(i);
      xmass5(i) = xmfrz(i);
    }
  }
  // Clin-6/06/02-end
  // C
  FEM_DO_SAFE(i, 1, maxstr) {
    ataui(i) = 0e0;
    zt1(i) = 0e0;
    zt2(i) = 0e0;
    // Clin-4/25/03 add zt3(I) to track longitudinal positions of
    // partons/strings:
    zt3(i) = 0e0;
    np(i) = 0;
  }
  FEM_DO_SAFE(i, 1, mul) {
    istrg = lstrg1(i);
    // Clin-9/2015 to avoid Floating-Point Exception:
    // C         TAU7 = SQRT(FT5(I) ** 2 - GZ5(I) ** 2)
    diff2 = fem::pow2(ft5(i)) - fem::pow2(gz5(i));
    if (diff2 < 0e0) {
      write(6, star), "2:I,ft5,gz5,diff2=", i, ft5(i), gz5(i), diff2;
      tau7 = 1e-6;
    } else {
      tau7 = fem::dsqrt(diff2);
    }
    // C
    ataui(istrg) += tau7;
    zt1(istrg) += gx5(i);
    zt2(istrg) += gy5(i);
    zt3(istrg) += gz5(i);
    np(istrg)++;
  }
  // C
  nstr = nsp + nst + cmn.nsi;
  // C
  // Clin-7/03/01 correct averaging on transverse coordinates, no shift needed:
  if (isoft == 3 || isoft == 4 || isoft == 5) {
    FEM_DO_SAFE(i, 1, nstr) {
      if (np(i) != 0) {
        ataui(i) = ataui(i) / np(i);
        zt1(i) = zt1(i) / np(i);
        zt2(i) = zt2(i) / np(i);
        zt3(i) = zt3(i) / np(i);
      }
    }
    return;
  }
  // Clin-7/03/01-end
  // C
  FEM_DO_SAFE(i, 1, nstr) {
    if (np(i) != 0) {
      ataui(i) = ataui(i) / np(i);
      zt1(i) = zt1(i) / np(i);
      zt2(i) = zt2(i) / np(i);
      zt3(i) = zt3(i) / np(i);
    } else {
      if (i <= nsp) {
        j = i;
        zt1(i) = fem::dble(yp(1, j));
        zt2(i) = fem::dble(yp(2, j));
        zt3(i) = 0e0;
      } else if (i > nsp && i <= nsp + nst) {
        j = i - nsp;
        zt1(i) = fem::dble(yt(1, j));
        zt2(i) = fem::dble(yt(2, j));
        zt3(i) = 0e0;
      } else {
        j = i - nsp - nst;
        zt1(i) = 0.5e0 * fem::dble((yp(1, iasg(j, 1)) + yt(1, iasg(j, 2))));
        zt2(i) = 0.5e0 * fem::dble((yp(2, iasg(j, 1)) + yt(2, iasg(j, 2))));
        zt3(i) = 0e0;
      }
    }
  }
  // C
  // Cbz6/28/99 flow1
  bb = hint1(19);
  FEM_DO_SAFE(i, 1, nstr) {
    if (np(i) != 0) {
      shift = 0e0;
    } else {
      shift = 0.5e0 * fem::dble(bb);
    }
    if (i <= nsp) {
      zt1(i) += shift;
    } else if (i > nsp && i <= nsp + nst) {
      zt1(i) = zt1(i) - shift;
    }
  }
  // Cbz6/28/99 flow1 end
  // C
}

}  // namespace AMPT
