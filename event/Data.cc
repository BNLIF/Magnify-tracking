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
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TPolyMarker3D.h"
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
    T_rec = (TTree*)rootFile->Get("T_rec");
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

void Data::DrawDQDX()
{
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
    if (reduced_chi2->size()>0) {
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
    int index = data_cluster_map[cluster_id];

    TH2F *h_proj_u = (TH2F*)gROOT->FindObject("h_proj_u");
    TH2F *h_proj_v = (TH2F*)gROOT->FindObject("h_proj_v");
    TH2F *h_proj_w = (TH2F*)gROOT->FindObject("h_proj_w");

    TGraph *g_rec_u = (TGraph*)gROOT->FindObject("g_rec_u");
    TGraph *g_rec_v = (TGraph*)gROOT->FindObject("g_rec_v");
    TGraph *g_rec_w = (TGraph*)gROOT->FindObject("g_rec_w");

    TH2F *h_pred_u = (TH2F*)gROOT->FindObject("h_pred_u");
    TH2F *h_pred_v = (TH2F*)gROOT->FindObject("h_pred_v");
    TH2F *h_pred_w = (TH2F*)gROOT->FindObject("h_pred_w");

    if (h_proj_u) { delete h_proj_u;}
    if (h_proj_v) { delete h_proj_v;}
    if (h_proj_w) { delete h_proj_w;}

    if (g_rec_u) { delete g_rec_u;}
    if (g_rec_v) { delete g_rec_v;}
    if (g_rec_w) { delete g_rec_w;}

    if (h_pred_u) { delete h_pred_u;}
    if (h_pred_v) { delete h_pred_v;}
    if (h_pred_w) { delete h_pred_w;}

    h_proj_u = new TH2F("h_proj_u", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_proj_v = new TH2F("h_proj_v", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_proj_w = new TH2F("h_proj_w", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);

    h_pred_u = new TH2F("h_pred_u", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
    h_pred_v = new TH2F("h_pred_v", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);
    h_pred_w = new TH2F("h_pred_w", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);

    TH2F *h[3] = {h_proj_u, h_proj_v, h_proj_w};
    TH2F *hpred[3] = {h_pred_u, h_pred_v, h_pred_w};
    int size = data_channel->at(index).size();
    TH2F *hc = 0;
    TH2F *hp = 0;
    for (int i=0; i<size; i++) {
        int x = data_channel->at(index).at(i);
        int y = data_time_slice->at(index).at(i);
        int z = data_charge->at(index).at(i);
        double charge_pred = data_charge_pred->at(index).at(i);
        double z_pred = (charge_pred-z+0.01)/(z+0.01);
        if (sign == 0) {
            z_pred = TMath::Abs(z_pred);
        }
        else if (sign>0) {
            if (z_pred<0) z_pred = 0;
        }
        else {
            if (z_pred>0) z_pred = 0;
            else z_pred = TMath::Abs(z_pred);
        }
        // double z_pred = TMath::Abs(charge_pred-z+0.01);
        // double z_pred = abs(charge_pred-z+0.01)/(data_charge_err->at(index).at(i)+0.01);
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

    int size2 = rec_u->at(currentCluster).size();
    g_rec_u = new TGraph(size2);
    g_rec_v = new TGraph(size2);
    g_rec_w = new TGraph(size2);
    TGraph *g[3] = {g_rec_u, g_rec_v, g_rec_w};

    for (int i=0; i<size2; i++) {
        double u = rec_u->at(currentCluster).at(i);
        double v = rec_v->at(currentCluster).at(i);
        double w = rec_w->at(currentCluster).at(i);
        double t = rec_t->at(currentCluster).at(i);
        // double q = rec_dQ->at(bin).at(i);
        g_rec_u->SetPoint(i, u, t);
        g_rec_v->SetPoint(i, v, t);
        g_rec_w->SetPoint(i, w, t);
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

    TH2F *h = (TH2F*)gROOT->FindObject("h_proj_u");
    if (zoomBin<0) {
        if (t0<-1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        }
        else {
            h->GetXaxis()->SetRangeUser(u0, u1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    }
    else {
        h->GetXaxis()->SetRangeUser(u-zoomBin, u+zoomBin);
        h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    }
    c1->GetPad(pad_proj)->Modified();
    c1->GetPad(pad_proj)->Update();

    h = (TH2F*)gROOT->FindObject("h_proj_v");
    if (zoomBin<0) {
        if (t0<-1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        }
        else {
            h->GetXaxis()->SetRangeUser(v0, v1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    }
    else {
        h->GetXaxis()->SetRangeUser(v-zoomBin, v+zoomBin);
        h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    }
    c1->GetPad(pad_proj+1)->Modified();
    c1->GetPad(pad_proj+1)->Update();

    h = (TH2F*)gROOT->FindObject("h_proj_w");
    if (zoomBin<0) {
        if (t0<-1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        }
        else {
            h->GetXaxis()->SetRangeUser(w0, w1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    }
    else {
        h->GetXaxis()->SetRangeUser(w-zoomBin, w+zoomBin);
        h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    }
    c1->GetPad(pad_proj+2)->Modified();
    c1->GetPad(pad_proj+2)->Update();

    // zoom prediction
    h = (TH2F*)gROOT->FindObject("h_pred_u");
    if (zoomBin<0) {
        if (t0<-1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        }
        else {
            h->GetXaxis()->SetRangeUser(u0, u1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    }
    else {
        h->GetXaxis()->SetRangeUser(u-zoomBin, u+zoomBin);
        h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    }
    c1->GetPad(pad_pred)->Modified();
    c1->GetPad(pad_pred)->Update();

    h = (TH2F*)gROOT->FindObject("h_pred_v");
    if (zoomBin<0) {
        if (t0<-1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        }
        else {
            h->GetXaxis()->SetRangeUser(v0, v1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    }
    else {
        h->GetXaxis()->SetRangeUser(v-zoomBin, v+zoomBin);
        h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    }
    c1->GetPad(pad_pred+1)->Modified();
    c1->GetPad(pad_pred+1)->Update();

    h = (TH2F*)gROOT->FindObject("h_pred_w");
    if (zoomBin<0) {
        if (t0<-1) {
            h->GetXaxis()->UnZoom();
            h->GetYaxis()->UnZoom();
        }
        else {
            h->GetXaxis()->SetRangeUser(w0, w1);
            h->GetYaxis()->SetRangeUser(t0, t1);
        }
    }
    else {
        h->GetXaxis()->SetRangeUser(w-zoomBin, w+zoomBin);
        h->GetYaxis()->SetRangeUser(t-zoomBin, t+zoomBin);
    }
    c1->GetPad(pad_pred+2)->Modified();
    c1->GetPad(pad_pred+2)->Update();

    // zoom 3D
    double x =  rec_x->at(currentCluster).at(pointIndex);
    double y =  rec_y->at(currentCluster).at(pointIndex);
    double z =  rec_z->at(currentCluster).at(pointIndex);
    TGraph2D *g = (TGraph2D*)gROOT->FindObject("g_3d");
    if (zoomBin<0) {
        g->GetXaxis()->UnZoom();
        g->GetYaxis()->UnZoom();
        g->GetZaxis()->UnZoom();
    }
    else {
        g->GetXaxis()->SetRangeUser(z-zoomBin*0.3, z+zoomBin*0.3);
        g->GetYaxis()->SetRangeUser(x-zoomBin*0.3, x+zoomBin*0.3);
        g->GetZaxis()->SetRangeUser(y-zoomBin*0.3, y+zoomBin*0.3);
    }
    c1->GetPad(pad_3d)->Modified();
    c1->GetPad(pad_3d)->Update();

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
    if (g_subclusters_3d.size()>0) {
        for (size_t i=0; i<g_subclusters_3d.size(); i++) {
            delete g_subclusters_3d.at(i);
        }
    }
    g_subclusters_3d.clear();


    int nSub = sub_id.size();
    cout << "#sub clusters: " << nSub << endl;
    for (int i=0; i<nSub; i++) {
        int nPoints = sub_end_index[i]-sub_start_index[i]+1;
        for (int k=0; k<3; k++) {
            pg[k]->push_back( new TGraph(nPoints) );
        }
        g_subclusters_3d.push_back( new TGraph2D(nPoints) );
        g_subclusters_3d[i]->SetName( TString::Format("g_3d_%i", sub_id[i]%1000) );

    }
    // cout << "g_subclusters_v: " << g_subclusters_v.size() << endl;


    int allPoints = rec_u->at(currentCluster).size();
    int i = 0;
    int currentSub = 0;
    int currentPointInSub = 0;
    while (i < allPoints) {
        double u = rec_u->at(currentCluster).at(i);
        double v = rec_v->at(currentCluster).at(i);
        double w = rec_w->at(currentCluster).at(i);
        double t = rec_t->at(currentCluster).at(i);
        int id = sub_cluster_id->at(currentCluster).at(i);

        g_subclusters_u[currentSub]->SetPoint(currentPointInSub, u, t);
        g_subclusters_v[currentSub]->SetPoint(currentPointInSub, v, t);
        g_subclusters_w[currentSub]->SetPoint(currentPointInSub, w, t);

        g_subclusters_3d[currentSub]->SetPoint(currentPointInSub,
            rec_z->at(currentCluster).at(i),
            rec_x->at(currentCluster).at(i),
            rec_y->at(currentCluster).at(i)
        );

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
        g_subclusters_3d[i]->SetLineWidth(2);
        g_subclusters_3d[i]->SetMarkerSize(0.5);
        g_subclusters_3d[i]->SetLineColor(colors[i%NC]);
        g_subclusters_3d[i]->SetMarkerStyle(markers[i%NM]);
        g_subclusters_3d[i]->Draw("pLINE,same");
    }  
    
    TLegend *leg = new TLegend(0.15, 0.50, 0.87, 0.87);
    for (int i=1; i<nSub; i++) {
        leg->AddEntry(pg[1]->at(i), TString::Format(" %i", sub_id[i]%1000), "lp");
    }
    c1->cd(pad_dqdx+1);    
    TH2F *hInfo = (TH2F*)gROOT->FindObject("hInfo");
    if (hInfo) { delete hInfo;}
    hInfo = new TH2F("hInfo","Info", 10, 0, 1, 10, 0, 1);
    hInfo->Draw();
    leg->SetNColumns(5);
    leg->Draw();

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
    TGraph2D *g = (TGraph2D*)gROOT->FindObject("g_3d");
    if (g) {delete g;}

    TGraph2D *gt = (TGraph2D*)gROOT->FindObject("gt_3d");
    if (gt) {delete gt;}

    int nPoints = rec_x->at(currentCluster).size();
    g = new TGraph2D(nPoints);
    g->SetName("g_3d");
    g->SetTitle("");
    for (int i=0; i<nPoints; i++) {
        g->SetPoint(i,
            rec_z->at(currentCluster).at(i),
            rec_x->at(currentCluster).at(i),
            rec_y->at(currentCluster).at(i)
        );
    }
    g->GetXaxis()->SetTitle("z [cm]");
    g->GetYaxis()->SetTitle("x [cm]");
    g->GetZaxis()->SetTitle("y [cm]");



    if (!isData){
      nPoints = true_x->at(0).size();
      gt = new TGraph2D(nPoints);
      gt->SetName("gt_3d");
      gt->SetTitle("");
      for (int i=0; i<nPoints; i++) {
        gt->SetPoint(i,
            true_z->at(0).at(i),
            true_x->at(0).at(i),
            true_y->at(0).at(i)
        );
      }
      gt->SetLineColor(2);
      gt->GetXaxis()->SetTitle("z [cm]");
      gt->GetYaxis()->SetTitle("x [cm]");
      gt->GetZaxis()->SetTitle("y [cm]");
    }


    c1->cd(pad_3d);
    // g->Draw("pLINE");
    g->SetMarkerSize(0);
    g->Draw("p");

    if(!isData){
        gt->Draw("LINEsame");
        TPolyMarker3D *gm = new TPolyMarker3D();
        gm->SetPoint(0,rec_z->at(currentCluster).front(),rec_x->at(currentCluster).front(),rec_y->at(currentCluster).front());
        gm->SetMarkerColor(4);
        gm->SetMarkerStyle(29);
        gm->SetMarkerSize(2);
        gm->Draw("*same");
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
            if (TMath::Abs(x-xx)<20 && TMath::Abs(y-yy)<20) {
                foundClusterIndex = i;
            }
        }
    }
    return foundClusterIndex;

}

void Data::DrawProjAll(int t0, int t1, int u0, int u1, int v0, int v1, int w0, int w1)
{
    int cluster_id = rec_cluster_id->at(currentCluster);
    int index = data_cluster_map[cluster_id];

    TH2F *h_proj_u_all = (TH2F*)gROOT->FindObject("h_proj_u_all");
    TH2F *h_proj_v_all = (TH2F*)gROOT->FindObject("h_proj_v_all");
    TH2F *h_proj_w_all = (TH2F*)gROOT->FindObject("h_proj_w_all");

    TGraph *g_rec_u_all = (TGraph*)gROOT->FindObject("g_rec_u_all");
    TGraph *g_rec_v_all = (TGraph*)gROOT->FindObject("g_rec_v_all");
    TGraph *g_rec_w_all = (TGraph*)gROOT->FindObject("g_rec_w_all");

    TH2F *h_pred_u_all = (TH2F*)gROOT->FindObject("h_pred_u_all");
    TH2F *h_pred_v_all = (TH2F*)gROOT->FindObject("h_pred_v_all");
    TH2F *h_pred_w_all = (TH2F*)gROOT->FindObject("h_pred_w_all");

    if (!h_proj_u_all) { h_proj_u_all = new TH2F("h_proj_u_all", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime); }
    if (!h_proj_v_all) { h_proj_v_all = new TH2F("h_proj_v_all", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime); }
    if (!h_proj_w_all) { h_proj_w_all = new TH2F("h_proj_w_all", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);}

    if (!h_pred_u_all) { h_pred_u_all = new TH2F("h_pred_u_all", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);}
    if (!h_pred_v_all) { h_pred_v_all = new TH2F("h_pred_v_all", "", nChannel_v, nChannel_u, nChannel_u+nChannel_v, nTime, 0, nTime);}
    if (!h_pred_w_all) { h_pred_w_all = new TH2F("h_pred_w_all", "", nChannel_w, nChannel_u+nChannel_v, nChannel_u+nChannel_v+nChannel_w, nTime, 0, nTime);}
    

    TH2F *h[3] = {h_proj_u_all, h_proj_v_all, h_proj_w_all};
    TH2F *hpred[3] = {h_pred_u_all, h_pred_v_all, h_pred_w_all};
    // 
    int nCluster = data_channel->size();
    TH2F *hc = 0;
    TH2F *hp = 0;
    for (int i=0; i<nCluster; i++) {
        int size = data_channel->at(i).size();
        for (int j=0; j<size; j++) {
            int x = data_channel->at(i).at(j);
            int y = data_time_slice->at(i).at(j);
            int z = data_charge->at(i).at(j);
            double charge_pred = data_charge_pred->at(i).at(j);
            double z_pred = (charge_pred-z+0.01)/(z+0.01);
            if (sign == 0) {
                z_pred = TMath::Abs(z_pred);
            }
            else if (sign>0) {
                if (z_pred<0) z_pred = 0;
            }
            else {
                if (z_pred>0) z_pred = 0;
                else z_pred = TMath::Abs(z_pred);
            }
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
    cout << "nCluster: " << nCluster << endl;
    cout << "nRec: " << nRec << endl;
    cout << "nPoint: " << nPoint << endl;
    if (!g_rec_u_all) { g_rec_u_all = new TGraph(nPoint);}
    if (!g_rec_v_all) { g_rec_v_all = new TGraph(nPoint);}
    if (!g_rec_w_all) { g_rec_w_all = new TGraph(nPoint);}

    TGraph *g[3] = {g_rec_u_all, g_rec_v_all, g_rec_w_all};
    int thisPoint = 0;
    for (int i=0; i<nRec; i++) {
        int size = rec_u->at(i).size();
        for (int j=0; j<size; j++) {
            double u = rec_u->at(i).at(j);
            double v = rec_v->at(i).at(j);
            double w = rec_w->at(i).at(j);
            double t = rec_t->at(i).at(j);
            g_rec_u_all->SetPoint(thisPoint, u, t);
            g_rec_v_all->SetPoint(thisPoint, v, t);
            g_rec_w_all->SetPoint(thisPoint, w, t);
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
}
