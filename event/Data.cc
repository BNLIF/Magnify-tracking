#include "Data.h"

#include "TH2F.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TColor.h"
#include "TBox.h"
#include "TPad.h"
#include "TPaletteAxis.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TLatex.h"


#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

Data::Data()
{}

double Data::applySign(double residual, int sign)
{
    if (sign == 0) {
        return TMath::Abs(residual);
    }
    if (sign > 0) {
        return residual < 0 ? 0.0 : residual;
    }
    // sign < 0
    return residual > 0 ? 0.0 : TMath::Abs(residual);
}

void Data::EnsureProjHistos()
{
    if (h_proj_u) return;
    h_proj_u = new TH2F("h_proj_u", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_proj_v = new TH2F("h_proj_v", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_proj_w = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);
    h_pred_u = new TH2F("h_pred_u", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_pred_v = new TH2F("h_pred_v", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_pred_w = new TH2F("h_pred_w", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);
}

void Data::EnsureProjAllHistos()
{
    if (h_proj_u_all) return;
    h_proj_u_all = new TH2F("h_proj_u_all", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_proj_v_all = new TH2F("h_proj_v_all", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_proj_w_all = new TH2F("h_proj_w_all", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);
    h_pred_u_all = new TH2F("h_pred_u_all", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_pred_v_all = new TH2F("h_pred_v_all", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_pred_w_all = new TH2F("h_pred_w_all", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);
}

Data::Data(const char* filename, int sign)
{
    this->sign = sign;
    c1 = 0;
    pad_dqdx = 1;
    pad_mc_compare = 2;
    pad_proj = 4;
    pad_3d = 3;
    pad_pred = 7;
    currentCluster = 0;
    doDrawBadCh = false;
    doDrawTrack = true;
    infoText = new TLatex(0.08, 0.20, "");

    rootFile = 0;
    T_true = 0;
    T_rec = 0;
    T_proj = 0;
    T_proj_data = 0;
    T_bad_ch = 0;

    nChannel_u = 2400;
    nChannel_v = 2400;
    nChannel_w = 3456;
    nTime = 2400;

    rec_cluster_id = new vector<int>;
    rec_x  = new vector<vector<double> >;
    rec_y  = new vector<vector<double> >;
    rec_z  = new vector<vector<double> >;
    rec_dQ = new vector<vector<double> >;
    rec_dx = new vector<vector<double> >;
    rec_L  = new vector<vector<double> >;
    rec_u  = new vector<vector<double> >;
    rec_v  = new vector<vector<double> >;
    rec_w  = new vector<vector<double> >;
    rec_t  = new vector<vector<double> >;
    reduced_chi2 = new vector<vector<double> >;
    com_dis  = new vector<vector<double> >;
    com_dtheta  = new vector<vector<double> >;
    true_dQ  = new vector<vector<double> >;
    flag_vertex  = new vector<vector<int> >;
    sub_cluster_id  = new vector<vector<int> >;

    data_cluster_id = new vector<int>;
    data_channel = new vector<vector<int> >;
    data_time_slice = new vector<vector<int> >;
    data_charge = new vector<vector<int> >;
    data_charge_pred = new vector<vector<int> >;;
    data_charge_err = new vector<vector<int> >;;

    true_x = new vector<vector<double> > ;
    true_y = new vector<vector<double> > ;
    true_z = new vector<vector<double> > ;

    stat_beg_dis = new vector<double>;
    stat_end_dis = new vector<double>;

    for (int i=0; i<3; i++) {
        currentPoint[i] = new TMarker(0, 0, 24);
        currentPoint[i]->SetMarkerColor(6);
        currentPoint[i]->SetMarkerSize(4);
    }
    dqdxPoint = new TMarker(0, 0, 24);
    dqdxPoint->SetMarkerColor(6);
    dqdxPoint->SetMarkerSize(4);

    h_3d_frame  = nullptr;
    g_3d_points = nullptr;
    gt_3d_line  = nullptr;
    g_3d_start  = nullptr;

    h_proj_u = nullptr; h_proj_v = nullptr; h_proj_w = nullptr;
    h_pred_u = nullptr; h_pred_v = nullptr; h_pred_w = nullptr;
    h_proj_u_all = nullptr; h_proj_v_all = nullptr; h_proj_w_all = nullptr;
    h_pred_u_all = nullptr; h_pred_v_all = nullptr; h_pred_w_all = nullptr;

    currentPoint3d = new TPolyMarker3D(1);
    currentPoint3d->SetMarkerStyle(24);
    currentPoint3d->SetMarkerColor(6);
    currentPoint3d->SetMarkerSize(2);

    LoadData(filename);
}



void Data::LoadData(const char* filename)
{
    rootFile = TFile::Open(filename);
    if (!rootFile) {
        string msg = "Unable to open "; msg += filename;
        throw runtime_error(msg.c_str());
    }

    T_true = (TTree*)rootFile->Get("T_true");
    T_rec = (TTree*)rootFile->Get("T_rec"); // modifyed from T_rec
    T_proj = (TTree*)rootFile->Get("T_proj");
    T_proj_data = (TTree*)rootFile->Get("T_proj_data");
    T_bad_ch = (TTree*)rootFile->Get("T_bad_ch");

    isData = T_true? false : true;
    cout << "loading data? " << isData << endl;

    LoadRec();
    LoadProj();
    LoadBadCh();
    if (!isData)
      LoadTruth();

}

void Data::LoadTruth(){
  T_true->SetBranchAddress("x",&true_x);
  T_true->SetBranchAddress("y",&true_y);
  T_true->SetBranchAddress("z",&true_z);
  T_true->GetEntry(0);
}

void Data::LoadRec()
{
    T_rec->SetBranchAddress("rec_x", &rec_x);
    T_rec->SetBranchAddress("rec_y", &rec_y);
    T_rec->SetBranchAddress("rec_z", &rec_z);
    T_rec->SetBranchAddress("rec_dQ", &rec_dQ);
    T_rec->SetBranchAddress("rec_dx", &rec_dx);
    T_rec->SetBranchAddress("rec_L", &rec_L);
    T_rec->SetBranchAddress("rec_cluster_id", &rec_cluster_id);
    T_rec->SetBranchAddress("rec_u", &rec_u);
    T_rec->SetBranchAddress("rec_v", &rec_v);
    T_rec->SetBranchAddress("rec_w", &rec_w);
    T_rec->SetBranchAddress("rec_t", &rec_t);
    if (T_rec->GetBranch("reduced_chi2")) {
        T_rec->SetBranchAddress("reduced_chi2", &reduced_chi2);
    }
    if (!isData) {
        T_rec->SetBranchAddress("com_dis", &com_dis);
        T_rec->SetBranchAddress("com_dtheta", &com_dtheta);
        T_rec->SetBranchAddress("true_dQ", &true_dQ);

	    T_rec->SetBranchAddress("stat_beg_dis",&stat_beg_dis);
	    T_rec->SetBranchAddress("stat_end_dis",&stat_end_dis);
    }
    T_rec->SetBranchAddress("flag_vertex", &flag_vertex);
    T_rec->SetBranchAddress("sub_cluster_id", &sub_cluster_id);

    T_rec->GetEntry(0);
    nCluster = rec_cluster_id->size();
    cout << " rec cluster id: ";
    for (int i=0; i<nCluster; i++) {
        cout << rec_cluster_id->at(i) << " ";
        int id = rec_cluster_id->at(i);
        rec_cluster_map[id] = i;
    }
    cout << endl;

    if (!isData){
      cout << "Starting point displacement: " << stat_beg_dis->at(0) << " cm" << std::endl;
      cout << "Ending point displacement  : " << stat_end_dis->at(0) << " cm" << std::endl;
    }

}

void Data::LoadProj()
{
    if (!T_proj_data) {
        cout << "no T_proj_data tree in the file ..." << endl;
        return;
    }
    T_proj_data->SetBranchAddress("cluster_id", &data_cluster_id);
    T_proj_data->SetBranchAddress("channel", &data_channel);
    T_proj_data->SetBranchAddress("time_slice", &data_time_slice);
    T_proj_data->SetBranchAddress("charge", &data_charge);
    T_proj_data->SetBranchAddress("charge_pred", &data_charge_pred);
    T_proj_data->SetBranchAddress("charge_err", &data_charge_err);

    T_proj_data->GetEntry(0);
    int size = data_cluster_id->size();
    cout << "data cluster id: ";
    for (int i=0; i<size; i++) {
        int id = data_cluster_id->at(i);
        data_cluster_map[id] = i;
        cout << id << " ";
    }
    cout << endl;
}

void Data::LoadBadCh()
{
    if (!T_bad_ch) {
        cout << "no bad channles in the file ..." << endl;
        return;
    }
    int chid, start, end;
    TLine *line = 0;
    T_bad_ch->SetBranchAddress("chid", &chid);
    T_bad_ch->SetBranchAddress("start_time", &start);
    T_bad_ch->SetBranchAddress("end_time", &end);
    int nEntries = T_bad_ch->GetEntries();
    double rebin = 4;
    for (int i=0; i<nEntries; i++) {
        T_bad_ch->GetEntry(i);
        bad_id.push_back(chid);
        bad_start.push_back(start);
        bad_end.push_back(end);
        line = new TLine(chid, start/rebin, chid, end/rebin);
        line->SetLineColorAlpha(kGray, 0.5);
        bad_lines.push_back(line);
    }
}

void Data::Clear()
{
    sub_id.clear();
    sub_start_index.clear();
    sub_end_index.clear();
}

void Data::DrawDQDX()
{
    Clear();
    
    TGraph *g = (TGraph*)gROOT->FindObject("g_dqdx");
    if (g) {
        delete g;
    }
    TVirtualPad *pad = c1->GetPad(pad_dqdx);

    int size = rec_dQ->at(currentCluster).size();
    sub_id.push_back(-1);
    sub_start_index.push_back(0);
    sub_end_index.push_back(0);
    int current_sub_index = 0;
    for (int i=0; i<size; i++) {
        int id = sub_cluster_id->at(currentCluster).at(i);
        if (id == sub_id[current_sub_index]) {
            sub_end_index[current_sub_index] = i;
        }
        else {
            sub_id.push_back(id);
            sub_start_index.push_back(i);
            sub_end_index.push_back(i);
            current_sub_index++;
        }
    }

    g = new TGraph(size);    
    for (int i=0; i<size; i++) {
      // std::cout << i << " " << rec_L->at(currentCluster).at(i) << std::endl;
        g->SetPoint(i, rec_L->at(currentCluster).at(i),
            rec_dQ->at(currentCluster).at(i)/1000/rec_dx->at(currentCluster).at(i));
    }
    g->SetName("g_dqdx");
    g->SetTitle(TString::Format("cluster id %i (index %i)", rec_cluster_id->at(currentCluster), currentCluster));
    g->GetXaxis()->SetTitle("Distance from start  [cm]");
    g->GetYaxis()->SetTitle("dQ/dx [1000 e^{-}/cm]");
    g->GetYaxis()->SetRangeUser(0,250);
    g->SetEditable(false);
    c1->cd(pad_dqdx);
    g->Draw("ALP");

    for (size_t i=0; i<sub_id.size(); i++) {
        cout << sub_id[i] << ", " << sub_start_index[i] << ", " << sub_end_index[i] << endl;
        double xloc = g->GetX()[sub_end_index[i]];
        TLine *line = new TLine(xloc+0.2, 0, xloc+0.5, 250);
        line->SetLineColorAlpha(kRed, 0.5);
        line->Draw();

        TLatex *tex = new TLatex(xloc, 255, TString::Format("%d", sub_id[i]%1000));
        tex->SetTextColor(kRed);
        tex->Draw();
    }
    

    if (!isData) {
        TGraph *g_dqdx_true = (TGraph*)gROOT->FindObject("g_dqdx_true");
        if (g_dqdx_true) {
            delete g_dqdx_true;
        }
        g_dqdx_true = new TGraph(size);

        for (int i=0; i<size; i++) {
            g_dqdx_true->SetPoint(i, rec_L->at(currentCluster).at(i),
                true_dQ->at(currentCluster).at(i)/1000/rec_dx->at(currentCluster).at(i));
        }
        g_dqdx_true->SetName("g_dqdx_true");
        g_dqdx_true->SetLineColor(kRed);
        g_dqdx_true->Draw("Lsame");
    }

    TGraph *g_reduced_chi2 = (TGraph*)gROOT->FindObject("g_reduced_chi2");
    if (reduced_chi2->size() > (size_t)currentCluster) {
        if (g_reduced_chi2) {
            delete g_reduced_chi2;
        }
        size = reduced_chi2->at(currentCluster).size();
        g_reduced_chi2 = new TGraph(size);

        for (int i=0; i<size; i++) {
            g_reduced_chi2->SetPoint(i, rec_L->at(currentCluster).at(i), reduced_chi2->at(currentCluster).at(i)*10);
        }
        g_reduced_chi2->SetName("g_reduced_chi2");
        g_reduced_chi2->SetTitle(TString::Format("cluster %i", rec_cluster_id->at(currentCluster)));
        g_reduced_chi2->SetLineColor(kBlue);
        g_reduced_chi2->SetMarkerColor(kBlue);
        g_reduced_chi2->SetMarkerStyle(7);
        // g_reduced_chi2->GetYaxis()->SetTitle("Distance (Data - MC) [cm]");
        g_reduced_chi2->Draw("LPsame");
        // // TGaxis *axis = new TGaxis(pad->GetUxmax(), pad->GetUymin(), pad->GetUxmax(), pad->GetUymax(),
        // //     pad->GetUymin(), pad->GetUymax(), 510,"+L");
        // // axis->SetLineColor(kRed);
        // // axis->Draw();
        TLegend *leg = new TLegend(0.67, 0.67, 0.87, 0.87);
        TLegendEntry *l1 = leg->AddEntry(g, "dQ/dx", "l");
        l1->SetTextColor(kBlack);
        TLegendEntry *l2 = leg->AddEntry(g_reduced_chi2, "#chi_{red.} x10", "l");
        l2->SetTextColor(kBlue);
        // leg->Draw();
    }


    pad->SetGridx();
    pad->SetGridy();
    pad->Modified();
    pad->Update();
}

void Data::DrawMCCompare()
{
    if (isData) return;
    c1->cd(pad_mc_compare);
    TVirtualPad *pad = c1->GetPad(pad_mc_compare);

    TGraph *g_com_dtheta = (TGraph*)gROOT->FindObject("g_com_dtheta");
    if (g_com_dtheta) {
        delete g_com_dtheta;
    }
    int size = com_dtheta->at(currentCluster).size();
    g_com_dtheta = new TGraph(size);

    for (int i=0; i<size; i++) {
        g_com_dtheta->SetPoint(i, rec_L->at(currentCluster).at(i), com_dtheta->at(currentCluster).at(i));
    }
    g_com_dtheta->SetName("g_com_dtheta");
    g_com_dtheta->SetTitle("Difference to MC");
    g_com_dtheta->GetXaxis()->SetTitle("Distance from start [cm]");
    g_com_dtheta->GetYaxis()->SetTitle("");
    g_com_dtheta->GetYaxis()->SetRangeUser(0,0.6);
    g_com_dtheta->Draw("ALP");

    TGraph *g = (TGraph*)gROOT->FindObject("g_com_dis");
    if (g) {
        delete g;
    }
    size = com_dis->at(currentCluster).size();
    g = new TGraph(size);

    for (int i=0; i<size; i++) {
        g->SetPoint(i, rec_L->at(currentCluster).at(i), com_dis->at(currentCluster).at(i));
    }
    g->SetName("g_com_dis");
    g->SetTitle(TString::Format("cluster %i", rec_cluster_id->at(currentCluster)));
    g->SetLineColor(kRed);
    g->SetMarkerColor(kRed);
    // g->GetYaxis()->SetTitle("Distance (Data - MC) [cm]");
    g->Draw("LPsame");
    // TGaxis *axis = new TGaxis(pad->GetUxmax(), pad->GetUymin(), pad->GetUxmax(), pad->GetUymax(),
    //     pad->GetUymin(), pad->GetUymax(), 510,"+L");
    // axis->SetLineColor(kRed);
    // axis->Draw();
    TLegend *leg = new TLegend(0.67, 0.67, 0.87, 0.87);
    TLegendEntry *l1 = leg->AddEntry(g_com_dtheta, "angle", "l");
    l1->SetTextColor(kBlack);
    TLegendEntry *l2 = leg->AddEntry(g, "distance", "l");
    l2->SetTextColor(kRed);
    leg->Draw();

    pad->SetGridx();
    pad->SetGridy();
    pad->Modified();
    pad->Update();

}


void Data::DrawProj()
{
    int cluster_id = rec_cluster_id->at(currentCluster);
    auto it = data_cluster_map.find(cluster_id);
    if (it == data_cluster_map.end()) {
        cout << "DrawProj: no projection data for cluster id " << cluster_id << endl;
        return;
    }
    int index = it->second;

    EnsureProjHistos();
    h_proj_u->Reset("ICES");
    h_proj_v->Reset("ICES");
    h_proj_w->Reset("ICES");
    h_pred_u->Reset("ICES");
    h_pred_v->Reset("ICES");
    h_pred_w->Reset("ICES");

    TGraph *g_rec_u = (TGraph*)gROOT->FindObject("g_rec_u");
    TGraph *g_rec_v = (TGraph*)gROOT->FindObject("g_rec_v");
    TGraph *g_rec_w = (TGraph*)gROOT->FindObject("g_rec_w");
    if (g_rec_u) { delete g_rec_u;}
    if (g_rec_v) { delete g_rec_v;}
    if (g_rec_w) { delete g_rec_w;}

    TH2F *h[3] = {h_proj_u, h_proj_v, h_proj_w};
    TH2F *hpred[3] = {h_pred_u, h_pred_v, h_pred_w};
    const auto& chs   = data_channel->at(index);
    const auto& ts    = data_time_slice->at(index);
    const auto& zs    = data_charge->at(index);
    const auto& zpres = data_charge_pred->at(index);
    int size = chs.size();
    TH2F *hc = 0;
    TH2F *hp = 0;
    for (int i=0; i<size; i++) {
        int x = chs[i];
        int y = ts[i];
        int z = zs[i];
        double charge_pred = zpres[i];
        double z_pred = applySign((charge_pred-z+0.01)/(z+0.01), sign);
        if (x<nChannel_u) {
            hc = h_proj_u;
            hp = h_pred_u;
        }
        else if (x<nChannel_u+nChannel_v) {
            hc = h_proj_v;
            hp = h_pred_v;
            x -= nChannel_u;
        }
        else {
            hc = h_proj_w;
            hp = h_pred_w;
            x -= (nChannel_u+nChannel_v);
        }
        hc->SetBinContent(x+1, y+1, z);
        hp->SetBinContent(x+1, y+1, z_pred);
    }

    const auto& us = rec_u->at(currentCluster);
    const auto& vs = rec_v->at(currentCluster);
    const auto& ws = rec_w->at(currentCluster);
    const auto& tts = rec_t->at(currentCluster);
    int size2 = us.size();
    g_rec_u = new TGraph(size2);
    g_rec_v = new TGraph(size2);
    g_rec_w = new TGraph(size2);
    TGraph *g[3] = {g_rec_u, g_rec_v, g_rec_w};

    for (int i=0; i<size2; i++) {
        double t = tts[i];
        g_rec_u->SetPoint(i, us[i], t);
        g_rec_v->SetPoint(i, vs[i], t);
        g_rec_w->SetPoint(i, ws[i], t);
    }
    for (int i=0; i<3; i++) {
        h[i]->SetTitle("Measured Charge");
        h[i]->GetXaxis()->SetTitle("Channel");
        h[i]->GetYaxis()->SetTitle("Time Slice");
        h[i]->GetZaxis()->SetRangeUser(500, 20000);
        int pad = pad_proj+i;
        c1->cd(pad);
        h[i]->Draw("colz");
        g[i]->SetLineWidth(2);
        g[i]->SetLineColor(6);
        if (doDrawTrack) {
            // g[i]->Draw("LPsame");
            // DrawSubclusters();
        }
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();

        hpred[i]->SetTitle("abs(charge_pred-charge)/charge");
        hpred[i]->GetXaxis()->SetTitle("Channel");
        hpred[i]->GetYaxis()->SetTitle("Time Slice");
        hpred[i]->GetZaxis()->SetRangeUser(0.01, 1);
        // hpred[i]->GetZaxis()->SetRangeUser(-1, 1);
        // hpred[i]->GetZaxis()->SetRangeUser(100, 2000);

        pad = pad_pred+i;
        c1->cd(pad);
        hpred[i]->Draw("colz");
        if (doDrawTrack) {
            // g[i]->Draw("LPsame");
            // DrawSubclusters();
        }
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();
    }

    if (doDrawTrack) {
        DrawSubclusters();
    }
    DrawBadCh();
}

static void zoomHisto(TH2F *h, int zoomBin, int t0, int t1, int x0, int x1, int xCenter, int tCenter)
{
    if (!h) return;
    if (zoomBin < 0) {
        if (t0 < -1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        } else {
            h->GetXaxis()->SetRangeUser(x0, x1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    } else {
        h->GetXaxis()->SetRangeUser(xCenter-zoomBin, xCenter+zoomBin);
        h->GetYaxis()->SetRangeUser(tCenter-zoomBin, tCenter+zoomBin);
    }
}

void Data::ZoomProj(int pointIndex, int zoomBin, int t0, int t1, int u0, int u1, int v0, int v1, int w0, int w1)
{
    int u =  rec_u->at(currentCluster).at(pointIndex);
    int v =  rec_v->at(currentCluster).at(pointIndex);
    int w =  rec_w->at(currentCluster).at(pointIndex);
    int t =  rec_t->at(currentCluster).at(pointIndex);

    cout << " index: " << pointIndex;
    cout << " u: " << u;
    cout << " v: " << v;
    cout << " w: " << w;
    cout << " t: " << t;
    cout << endl;

    zoomHisto(h_proj_u, zoomBin, t0, t1, u0, u1, u, t);
    c1->GetPad(pad_proj)->Modified();
    c1->GetPad(pad_proj)->Update();

    zoomHisto(h_proj_v, zoomBin, t0, t1, v0, v1, v, t);
    c1->GetPad(pad_proj+1)->Modified();
    c1->GetPad(pad_proj+1)->Update();

    zoomHisto(h_proj_w, zoomBin, t0, t1, w0, w1, w, t);
    c1->GetPad(pad_proj+2)->Modified();
    c1->GetPad(pad_proj+2)->Update();

    zoomHisto(h_pred_u, zoomBin, t0, t1, u0, u1, u, t);
    c1->GetPad(pad_pred)->Modified();
    c1->GetPad(pad_pred)->Update();

    zoomHisto(h_pred_v, zoomBin, t0, t1, v0, v1, v, t);
    c1->GetPad(pad_pred+1)->Modified();
    c1->GetPad(pad_pred+1)->Update();

    zoomHisto(h_pred_w, zoomBin, t0, t1, w0, w1, w, t);
    c1->GetPad(pad_pred+2)->Modified();
    c1->GetPad(pad_pred+2)->Update();

    // zoom 3D (ROOT-X=reco-z, ROOT-Y=reco-x, ROOT-Z=reco-y)
    {
        double xr = rec_x->at(currentCluster).at(pointIndex);
        double yr = rec_y->at(currentCluster).at(pointIndex);
        double zr = rec_z->at(currentCluster).at(pointIndex);
        if (h_3d_frame) {
            if (zoomBin<0) {
                h_3d_frame->GetXaxis()->UnZoom();
                h_3d_frame->GetYaxis()->UnZoom();
                h_3d_frame->GetZaxis()->UnZoom();
            } else {
                h_3d_frame->GetXaxis()->SetRangeUser(zr-zoomBin*0.3, zr+zoomBin*0.3);
                h_3d_frame->GetYaxis()->SetRangeUser(xr-zoomBin*0.3, xr+zoomBin*0.3);
                h_3d_frame->GetZaxis()->SetRangeUser(yr-zoomBin*0.3, yr+zoomBin*0.3);
            }
            c1->GetPad(pad_3d)->Modified();
            c1->GetPad(pad_3d)->Update();
        }
    }

    // zoom 1D
    TGraph *g_dqdx = (TGraph*)gROOT->FindObject("g_dqdx");
    if (zoomBin<0) {
        g_dqdx->GetXaxis()->UnZoom();
    }
    else {
        double xx, yy;
        g_dqdx->GetPoint(pointIndex, xx, yy);
        g_dqdx->GetXaxis()->SetRangeUser(xx-10, xx+10);
        c1->cd(pad_dqdx);
        dqdxPoint->SetX(xx);
        dqdxPoint->SetY(yy);
        dqdxPoint->Draw();
    }
    c1->GetPad(pad_dqdx)->Modified();
    c1->GetPad(pad_dqdx)->Update();

}

void Data::DrawSubclusters()
{
    vector<TGraph*>*  pg[3] = {&g_subclusters_u, &g_subclusters_v, &g_subclusters_w};
    for (int k=0; k<3; k++) {
        if (pg[k]->size()>0) {
            for (size_t i=0; i<pg[k]->size(); i++) {
                delete pg[k]->at(i);
            }
        }
        pg[k]->clear();
    }

    // 3D sub-cluster cleanup is handled by Draw3D (via pad->Clear("nodelete")).
    // Vectors are already empty when DrawSubclusters is called from DrawNewCluster.
    // On ToggleDrawTrack re-entry (DrawSubclusters called without Draw3D), the
    // vectors may be non-empty; clear them here using the same nodelete strategy.
    if (!g_subclusters_3d_lines.empty()) {
        c1->GetPad(pad_3d)->Clear("nodelete");
        for (size_t i=0; i<g_subclusters_3d_lines.size(); i++) delete g_subclusters_3d_lines[i];
        g_subclusters_3d_lines.clear();
        for (size_t i=0; i<g_subclusters_3d_marks.size(); i++) delete g_subclusters_3d_marks[i];
        g_subclusters_3d_marks.clear();
    }

    int nSub = sub_id.size();
    cout << "#sub clusters: " << nSub << endl;
    for (int i=0; i<nSub; i++) {
        int nPoints = sub_end_index[i]-sub_start_index[i]+1;
        for (int k=0; k<3; k++) {
            pg[k]->push_back( new TGraph(nPoints) );
        }
        g_subclusters_3d_lines.push_back(new TPolyLine3D(nPoints));
        g_subclusters_3d_marks.push_back(new TPolyMarker3D(nPoints));
    }

    int allPoints = rec_u->at(currentCluster).size();
    int i = 0;
    int currentSub = 0;
    int currentPointInSub = 0;
    while (i < allPoints && currentSub < nSub) {
        double u = rec_u->at(currentCluster).at(i);
        double v = rec_v->at(currentCluster).at(i);
        double w = rec_w->at(currentCluster).at(i);
        double t = rec_t->at(currentCluster).at(i);

        g_subclusters_u[currentSub]->SetPoint(currentPointInSub, u, t);
        g_subclusters_v[currentSub]->SetPoint(currentPointInSub, v, t);
        g_subclusters_w[currentSub]->SetPoint(currentPointInSub, w, t);

        double zr = rec_z->at(currentCluster).at(i);
        double xr = rec_x->at(currentCluster).at(i);
        double yr = rec_y->at(currentCluster).at(i);
        g_subclusters_3d_lines[currentSub]->SetPoint(currentPointInSub, zr, xr, yr);
        g_subclusters_3d_marks[currentSub]->SetPoint(currentPointInSub, zr, xr, yr);

        i++;
        currentPointInSub++;

        if (i>sub_end_index[currentSub]) {
            currentSub++;
            currentPointInSub = 0;
        }
    }
    const int NC = 7;
    int colors[] = {1, 6, 2, 4, 8, 7, 28};
    const int NM = 3;
    int markers[] = {32, 25, 26};

    for (int k=0; k<3; k++) {
        int pad = pad_proj+k;
        c1->cd(pad);
        for (int i=1; i<nSub; i++) {
            pg[k]->at(i)->SetLineWidth(2);
            pg[k]->at(i)->SetMarkerSize(0.5);
            pg[k]->at(i)->SetLineColor(colors[i%NC]);
            pg[k]->at(i)->SetMarkerStyle(markers[i%NM]);
            pg[k]->at(i)->Draw("LPsame");
        }
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();

        pad = pad_pred+k;
        c1->cd(pad);
        for (int i=1; i<nSub; i++) {
            pg[k]->at(i)->Draw("LPsame");
        }
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();
    }

    c1->cd(pad_3d);
    for (int i=1; i<nSub; i++) {
        g_subclusters_3d_lines[i]->SetLineWidth(2);
        g_subclusters_3d_lines[i]->SetLineColor(colors[i%NC]);
        g_subclusters_3d_lines[i]->Draw();

        g_subclusters_3d_marks[i]->SetMarkerSize(0.5);
        g_subclusters_3d_marks[i]->SetMarkerColor(colors[i%NC]);
        g_subclusters_3d_marks[i]->SetMarkerStyle(markers[i%NM]);
        g_subclusters_3d_marks[i]->Draw();
    }
    c1->GetPad(pad_3d)->Modified();
    c1->GetPad(pad_3d)->Update();

    if (isData) {
        // pad_mc_compare is unused for data files; use the full pad for
        // a sub-cluster id legend backed by an empty hInfo histogram.
        TLegend *leg = new TLegend(0.15, 0.40, 0.87, 0.87);
        for (int i=1; i<nSub; i++) {
            TLegendEntry* le = leg->AddEntry(pg[1]->at(i), TString::Format(" %i", sub_id[i]%1000), "p");
            le->SetTextColor(colors[i%NC]);
        }
        c1->cd(pad_mc_compare);
        TH2F *hInfo = (TH2F*)gROOT->FindObject("hInfo");
        if (hInfo) { delete hInfo;}
        hInfo = new TH2F("hInfo","Info", 10, 0, 1, 10, 0, 1);
        hInfo->Draw();
        leg->SetNColumns(5);
        leg->Draw();
        infoText->Draw();
    }
    else {
        // MC mode: pad_mc_compare shows the data-vs-truth comparison.
        // Overlay a compact sub-cluster id legend on pad_3d so we don't
        // clobber the MC compare panel.
        TLegend *leg = new TLegend(0.0, 0.78, 1.0, 1.0);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        for (int i=1; i<nSub; i++) {
            TLegendEntry* le = leg->AddEntry(pg[1]->at(i), TString::Format(" %i", sub_id[i]%1000), "p");
            le->SetTextColor(colors[i%NC]);
        }
        leg->SetNColumns(5);
        c1->cd(pad_3d);
        leg->Draw();
        c1->GetPad(pad_3d)->Modified();
        c1->GetPad(pad_3d)->Update();
    }

}

void Data::DrawPoint(int pointIndex)
{
    double x[3] = {0};
    x[0] =  rec_u->at(currentCluster).at(pointIndex);
    x[1] =  rec_v->at(currentCluster).at(pointIndex);
    x[2] =  rec_w->at(currentCluster).at(pointIndex);
    double t =  rec_t->at(currentCluster).at(pointIndex);

    for (int i=0; i<3; i++) {
        c1->cd(pad_proj+i);
        currentPoint[i]->SetX(x[i]);
        currentPoint[i]->SetY(t);
        currentPoint[i]->Draw();
        c1->GetPad(pad_proj+i)->Modified();
        c1->GetPad(pad_proj+i)->Update();

        c1->cd(pad_pred+i);
        currentPoint[i]->Draw();
        c1->GetPad(pad_proj+i)->Modified();
        c1->GetPad(pad_proj+i)->Update();
    }

    // highlight same point in 3D view
    if (currentPoint3d) {
        currentPoint3d->SetPoint(0,
            rec_z->at(currentCluster).at(pointIndex),
            rec_x->at(currentCluster).at(pointIndex),
            rec_y->at(currentCluster).at(pointIndex));
        c1->cd(pad_3d);
        currentPoint3d->Draw();
        c1->GetPad(pad_3d)->Modified();
        c1->GetPad(pad_3d)->Update();
    }

}

void Data::DrawBadCh()
{
    int size = bad_lines.size();
    for (int pad=pad_proj; pad<pad_proj+3; pad++) {
        c1->cd(pad);
        TVirtualPad *p = c1->GetPad(pad);
        for (int i=0; i<size; i++) {
            if ( (pad==pad_proj && bad_id[i]<2400)
                || (pad==pad_proj+1 && bad_id[i]>=2400 && bad_id[i]<4800)
                || (pad==pad_proj+2 && bad_id[i]>=4800)
            ) {
                if (doDrawBadCh) {
                    bad_lines[i]->Draw();
                }
                else {
                   p->GetListOfPrimitives()->Remove(bad_lines[i]);
                }

            }
        }
        p->Modified();
        p->Update();
    }
}

void Data::Draw3D()
{
    c1->cd(pad_3d);

    // Remove ALL primitives from pad_3d without deleting them. This preempts
    // TH3F::Paint's internal gPad->Clear() from double-freeing our objects
    // during the canvas repaint triggered by later DrawProj pad updates.
    c1->GetPad(pad_3d)->Clear("nodelete");

    // Safe to delete now: removed from pad by Clear above.
    if (h_3d_frame) { delete h_3d_frame; h_3d_frame = nullptr; }
    if (g_3d_points) { delete g_3d_points; g_3d_points = nullptr; }
    if (gt_3d_line)  { delete gt_3d_line;  gt_3d_line  = nullptr; }
    if (g_3d_start)  { delete g_3d_start;  g_3d_start  = nullptr; }

    // Sub-cluster 3D objects from previous DrawSubclusters call are also
    // removed from the pad by Clear above; delete them here.
    for (size_t i=0; i<g_subclusters_3d_lines.size(); i++) delete g_subclusters_3d_lines[i];
    g_subclusters_3d_lines.clear();
    for (size_t i=0; i<g_subclusters_3d_marks.size(); i++) delete g_subclusters_3d_marks[i];
    g_subclusters_3d_marks.clear();

    int nPoints = rec_x->at(currentCluster).size();

    // Compute bounding box. Coordinate mapping (preserved from original):
    //   ROOT-X = reco-z,  ROOT-Y = reco-x,  ROOT-Z = reco-y
    double xMin= 1e30, xMax=-1e30;
    double yMin= 1e30, yMax=-1e30;
    double zMin= 1e30, zMax=-1e30;
    for (int i=0; i<nPoints; i++) {
        double X = rec_z->at(currentCluster).at(i);
        double Y = rec_x->at(currentCluster).at(i);
        double Z = rec_y->at(currentCluster).at(i);
        if (X<xMin) xMin=X;  if (X>xMax) xMax=X;
        if (Y<yMin) yMin=Y;  if (Y>yMax) yMax=Y;
        if (Z<zMin) zMin=Z;  if (Z>zMax) zMax=Z;
    }
    if (!isData) {
        int nt = (int)true_x->at(0).size();
        for (int i=0; i<nt; i++) {
            double X = true_z->at(0).at(i);
            double Y = true_x->at(0).at(i);
            double Z = true_y->at(0).at(i);
            if (X<xMin) xMin=X;  if (X>xMax) xMax=X;
            if (Y<yMin) yMin=Y;  if (Y>yMax) yMax=Y;
            if (Z<zMin) zMin=Z;  if (Z>zMax) zMax=Z;
        }
    }
    // 5% margin on each axis; floor span to 1 cm so single-point cluster
    // still yields a visible frame.
    double frameMargin = 0.05;
    double dx = std::max(xMax-xMin, 1.0) * frameMargin;
    double dy = std::max(yMax-yMin, 1.0) * frameMargin;
    double dz = std::max(zMax-zMin, 1.0) * frameMargin;
    xMin -= dx; xMax += dx;
    yMin -= dy; yMax += dy;
    zMin -= dz; zMax += dz;

    // Empty TH3F frame: establishes TView3D on pad_3d and provides the axes
    // that ZoomProj calls SetRangeUser/UnZoom on. No bins filled -> no mesh.
    h_3d_frame = new TH3F("h_3d_frame", "",
                          10, xMin, xMax,
                          10, yMin, yMax,
                          10, zMin, zMax);
    h_3d_frame->SetDirectory(0);
    h_3d_frame->SetStats(0);
    h_3d_frame->GetXaxis()->SetTitle("z [cm]");
    h_3d_frame->GetYaxis()->SetTitle("x [cm]");
    h_3d_frame->GetZaxis()->SetTitle("y [cm]");
    h_3d_frame->Draw();

    // Base point cloud: small gray dots for all reco hits. The sub-cluster
    // overlay (drawn later by DrawSubclusters) provides the bright coloring.
    g_3d_points = new TPolyMarker3D(nPoints);
    for (int i=0; i<nPoints; i++) {
        g_3d_points->SetPoint(i,
            rec_z->at(currentCluster).at(i),
            rec_x->at(currentCluster).at(i),
            rec_y->at(currentCluster).at(i));
    }
    g_3d_points->SetMarkerStyle(1);
    g_3d_points->SetMarkerColor(15);
    g_3d_points->SetMarkerSize(0.3);
    g_3d_points->Draw();

    if (!isData) {
        int nt = (int)true_x->at(0).size();
        gt_3d_line = new TPolyLine3D(nt);
        for (int i=0; i<nt; i++) {
            gt_3d_line->SetPoint(i,
                true_z->at(0).at(i),
                true_x->at(0).at(i),
                true_y->at(0).at(i));
        }
        gt_3d_line->SetLineColor(2);
        gt_3d_line->SetLineWidth(2);
        gt_3d_line->Draw();

        g_3d_start = new TPolyMarker3D(1);
        g_3d_start->SetPoint(0,
            rec_z->at(currentCluster).front(),
            rec_x->at(currentCluster).front(),
            rec_y->at(currentCluster).front());
        g_3d_start->SetMarkerColor(4);
        g_3d_start->SetMarkerStyle(29);
        g_3d_start->SetMarkerSize(2);
        g_3d_start->Draw();
    }

    c1->GetPad(pad_3d)->Modified();
    c1->GetPad(pad_3d)->Update();
}

void Data::DrawNewCluster()
{
    DrawDQDX();
    DrawMCCompare();
    Draw3D();
    // DrawComDtheta();
    DrawProj();
    DrawPoint(0);
}

int Data::FindClusterIndex(double x, double y)
{
    int foundClusterIndex = -1;
    vector<vector<double> >* rec=0;
    if (x>0 && x<nChannel_u) { rec = rec_u; }
    else if (x>nChannel_u && x<nChannel_u+nChannel_v) { rec = rec_v; }
    else if (x>nChannel_u+nChannel_v && x<nChannel_u+nChannel_v+nChannel_w) { rec = rec_w; }

    if(!rec) { return -1; }

    int nRec = rec->size();
    for (int i=0; i<nRec; i++) {
        int size = rec->at(i).size();
        for (int j=0; j<size; j++) {
            double xx = rec->at(i).at(j);
            double yy = rec_t->at(i).at(j);
            if (TMath::Abs(x-xx)<10 && TMath::Abs(y-yy)<10) {
                foundClusterIndex = i;
            }
        }
    }
    return foundClusterIndex;

}

int Data::FindPointIndex(double x, double y)
{
    int foundIndex = -1;
    vector<vector<double> >* rec=0;
    if (x>0 && x<nChannel_u) { rec = rec_u; }
    else if (x>nChannel_u && x<nChannel_u+nChannel_v) { rec = rec_v; }
    else if (x>nChannel_u+nChannel_v && x<nChannel_u+nChannel_v+nChannel_w) { rec = rec_w; }

    if(!rec) { return -1; }

    const auto& xs = rec->at(currentCluster);
    const auto& ts = rec_t->at(currentCluster);
    int size = xs.size();
    double bestDist2 = 50.0; // 5x5 + 5x5 search window, squared
    for (int i=0; i<size; i++) {
        double dx = x - xs[i];
        double dy = y - ts[i];
        double d2 = dx*dx + dy*dy;
        if (d2 < bestDist2) {
            bestDist2 = d2;
            foundIndex = i;
        }
    }
    return foundIndex;

}

void Data::DrawProjAll(int t0, int t1, int u0, int u1, int v0, int v1, int w0, int w1)
{
    EnsureProjAllHistos();
    h_proj_u_all->Reset("ICES");
    h_proj_v_all->Reset("ICES");
    h_proj_w_all->Reset("ICES");
    h_pred_u_all->Reset("ICES");
    h_pred_v_all->Reset("ICES");
    h_pred_w_all->Reset("ICES");

    TGraph *g_rec_u_all = (TGraph*)gROOT->FindObject("g_rec_u_all");
    TGraph *g_rec_v_all = (TGraph*)gROOT->FindObject("g_rec_v_all");
    TGraph *g_rec_w_all = (TGraph*)gROOT->FindObject("g_rec_w_all");

    TH2F *h[3] = {h_proj_u_all, h_proj_v_all, h_proj_w_all};
    TH2F *hpred[3] = {h_pred_u_all, h_pred_v_all, h_pred_w_all};
    int nDataClusters = data_channel->size();
    TH2F *hc = 0;
    TH2F *hp = 0;
    for (int i=0; i<nDataClusters; i++) {
        const auto& chs   = data_channel->at(i);
        const auto& ts    = data_time_slice->at(i);
        const auto& zs    = data_charge->at(i);
        const auto& zpres = data_charge_pred->at(i);
        int size = chs.size();
        for (int j=0; j<size; j++) {
            int x = chs[j];
            int y = ts[j];
            int z = zs[j];
            double charge_pred = zpres[j];
            double z_pred = applySign((charge_pred-z+0.01)/(z+0.01), sign);
            if (x<nChannel_u) {
                hc = h_proj_u_all;
                hp = h_pred_u_all;
            }
            else if (x<nChannel_u+nChannel_v) {
                hc = h_proj_v_all;
                hp = h_pred_v_all;
                x -= nChannel_u;
            }
            else {
                hc = h_proj_w_all;
                hp = h_pred_w_all;
                x -= (nChannel_u+nChannel_v);
            }
            hc->SetBinContent(x+1, y+1, z);
            hp->SetBinContent(x+1, y+1, z_pred);
        }
    }

    int nPoint = 0;
    int nRec = rec_u->size();
    for (int i=0; i<nRec; i++) {
        nPoint += rec_u->at(i).size();
    }
    cout << "nDataClusters: " << nDataClusters << endl;
    cout << "nRec: " << nRec << endl;
    cout << "nPoint: " << nPoint << endl;
    if (!g_rec_u_all) { g_rec_u_all = new TGraph(nPoint);}
    if (!g_rec_v_all) { g_rec_v_all = new TGraph(nPoint);}
    if (!g_rec_w_all) { g_rec_w_all = new TGraph(nPoint);}

    TGraph *g[3] = {g_rec_u_all, g_rec_v_all, g_rec_w_all};
    int thisPoint = 0;
    for (int i=0; i<nRec; i++) {
        const auto& us = rec_u->at(i);
        const auto& vs = rec_v->at(i);
        const auto& ws = rec_w->at(i);
        const auto& tts = rec_t->at(i);
        int size = us.size();
        for (int j=0; j<size; j++) {
            double t = tts[j];
            g_rec_u_all->SetPoint(thisPoint, us[j], t);
            g_rec_v_all->SetPoint(thisPoint, vs[j], t);
            g_rec_w_all->SetPoint(thisPoint, ws[j], t);
            thisPoint++;
        }
    }
    // cout << "thisPoint: " << thisPoint << endl;

    int min[3] = {u0, v0, w0};
    int max[3] = {u1, v1, w1};
    
    for (int i=0; i<3; i++) {
        h[i]->SetTitle("Measured Charge");
        h[i]->GetXaxis()->SetTitle("Channel");
        h[i]->GetYaxis()->SetTitle("Time Slice");
        h[i]->GetZaxis()->SetRangeUser(500, 20000);

        if (t0<-1) {
            h[i]->GetXaxis()->UnZoom();
            h[i]->GetYaxis()->UnZoom();
        }
        else {
            h[i]->GetXaxis()->SetRangeUser(min[i], max[i]);
            h[i]->GetYaxis()->SetRangeUser(t0, t1);
        }

        int pad = pad_proj+i;
        c1->cd(pad);
        h[i]->Draw("colz");
        g[i]->SetLineWidth(2);
        g[i]->SetLineColor(6);
        if (doDrawTrack) {
            g[i]->Draw("Psame");
        }
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();

        hpred[i]->SetTitle("abs(charge_pred-charge)/charge");
        hpred[i]->GetXaxis()->SetTitle("Channel");
        hpred[i]->GetYaxis()->SetTitle("Time Slice");
        hpred[i]->GetZaxis()->SetRangeUser(0.01, 1);
        if (t0<-1) {
            hpred[i]->GetXaxis()->UnZoom();
            hpred[i]->GetYaxis()->UnZoom();
        }
        else {
            hpred[i]->GetXaxis()->SetRangeUser(min[i], max[i]);
            hpred[i]->GetYaxis()->SetRangeUser(t0, t1);
        }

        pad = pad_pred+i;
        c1->cd(pad);
        hpred[i]->Draw("colz");
        if (doDrawTrack) {
            g[i]->Draw("Psame");
        }
        c1->GetPad(pad)->Modified();
        c1->GetPad(pad)->Update();
    }

    DrawBadCh();
}




// void Data::Project()
// {

//     h_proj[0] = new TH2F("h_proj_u", "", nChannel_u, 1-0.5, nChannel_u-0.5, nTime, 1-0.5, nTime-0.5);
//     h_proj[1] = new TH2F("h_proj_v", "", nChannel_v, nChannel_u-0.5, nChannel_u+nChannel_v-0.5, nTime, 1-0.5, nTime-0.5);
//     h_proj[2] = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v-0.5, nChannel_u+nChannel_v+nChannel_w-0.5, nTime, 1-0.5, nTime-0.5);

//     T_proj_data->Project("h_proj_u", "time_slice:channel", "charge/500", "channel<2400");
//     T_proj_data->Project("h_proj_v", "time_slice:channel", "charge/500", "channel<2400 && channel<4800");
//     T_proj_data->Project("h_proj_w", "time_slice:channel", "charge/500", "channel>4800");

//     //creating boxes
//     int threshold = 0;
//     TBox *box = 0;
//     for (int plane=0; plane<3; plane++) {
//         for (int i=1; i<=nChannel_u; i++) {
//             for (int j=1; j<=nTime; j++) {
//                 double content = h_proj[plane]->GetBinContent(i, j);
//                 if (TMath::Abs(content)>threshold) {
//                     box = new TBox(
//                         h_proj[plane]->GetXaxis()->GetBinLowEdge(i),
//                         h_proj[plane]->GetYaxis()->GetBinLowEdge(j),
//                         h_proj[plane]->GetXaxis()->GetBinUpEdge(i),
//                         h_proj[plane]->GetYaxis()->GetBinUpEdge(j)
//                     );
//                     box->SetFillColor(kRed);
//                     boxes[plane].push_back(box);
//                     box_values[plane].push_back(TMath::Abs(content));
//                 }
//             }
//         }
//         cout << "plane " << plane << ": " << boxes[plane].size() <<  " boxes created. " << endl;
//     }

// }

// void Data::DrawProjection(int plane)
// {

//     h_proj[plane]->Clear();
//     h_proj[plane]->Draw("colz");
//     h_proj[plane]->GetZaxis()->SetRangeUser(10, 20);
//     pad->Update();

//     TPaletteAxis *palette = (TPaletteAxis*)h_proj[plane]->GetListOfFunctions()->FindObject("palette");
//     int size = boxes[plane].size();
//     for (int i=0; i<size; i++) {
//         boxes[plane][i]->SetFillColor(palette->GetValueColor(box_values[plane][i]));
//         boxes[plane][i]->Draw();
//     }
// }


Data::~Data()
{
    delete rootFile;
    delete h_3d_frame;
    delete g_3d_points;
    delete gt_3d_line;
    delete g_3d_start;
    delete currentPoint3d;
    delete infoText;
    delete dqdxPoint;
    delete h_proj_u; delete h_proj_v; delete h_proj_w;
    delete h_pred_u; delete h_pred_v; delete h_pred_w;
    delete h_proj_u_all; delete h_proj_v_all; delete h_proj_w_all;
    delete h_pred_u_all; delete h_pred_v_all; delete h_pred_w_all;
    for (int i=0; i<3; i++) delete currentPoint[i];
    for (size_t i=0; i<bad_lines.size(); i++) delete bad_lines[i];
    for (size_t i=0; i<g_subclusters_u.size(); i++) delete g_subclusters_u[i];
    for (size_t i=0; i<g_subclusters_v.size(); i++) delete g_subclusters_v[i];
    for (size_t i=0; i<g_subclusters_w.size(); i++) delete g_subclusters_w[i];
    for (size_t i=0; i<g_subclusters_3d_lines.size(); i++) delete g_subclusters_3d_lines[i];
    for (size_t i=0; i<g_subclusters_3d_marks.size(); i++) delete g_subclusters_3d_marks[i];
}
