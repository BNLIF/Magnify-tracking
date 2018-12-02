#include "GuiController.h"
#include "MainWindow.h"
#include "ViewWindow.h"
#include "ControlWindow.h"
#include "Data.h"

#include "TApplication.h"
#include "TSystem.h"
#include "TExec.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGFileDialog.h"

#include "TGMenu.h"
#include "TGNumberEntry.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TBox.h"
#include "TLine.h"
#include "TColor.h"

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


GuiController::GuiController(const TGWindow *p, int w, int h, const char* fn, int sign)
{
    mw = new MainWindow(p, w, h);
    vw = mw->fViewWindow;
    cw = mw->fControlWindow;
    currentCluster = 0;
    currentPointIndex = 0;

    TString filename;
    if (!fn) {
        filename = OpenDialog();
    }
    else {
        filename = fn;
    }
    data = new Data(filename.Data(), sign);
    data->c1 = vw->can;
    SetCurrentCluster(0);
    data->DrawNewCluster();


    // mw->SetWindowName(TString::Format("Magnify: run %i, sub-run %i, event %i",
    //     data->runNo, data->subRunNo, data->eventNo));

    // for (int i=0; i<6; i++) {
    //     vw->can->cd(i+1);
    //     data->wfs.at(i)->Draw2D();
    // }
    // for (int i=0; i<3; i++) {
    //     vw->can->cd(i+7);
    //     int chanNo = data->wfs.at(i)->firstChannel;
    //     data->wfs.at(i)->Draw1D(chanNo);
    //     TH1F *h = data->wfs.at(i+3)->Draw1D(chanNo, "same"); // draw calib
    //     h->SetLineColor(kRed);
    //     hCurrent[i] = h;
    // }



    InitConnections();

}

GuiController::~GuiController()
{
    // gApplication->Terminate(0);
}

void GuiController::InitConnections()
{
    mw->fMenuFile->Connect("Activated(int)", "GuiController", this, "HandleMenu(int)");

    cw->clusterEntry->SetLimitValues(0, data->nCluster-1);
    cw->clusterEntry->Connect("ValueSet(Long_t)", "GuiController", this, "ClusterChanged(int)");

    cw->zoomEntry->Connect("ValueSet(Long_t)", "GuiController", this, "ZoomChanged()");
    cw->badChanelButton->Connect("Clicked()", "GuiController", this, "ToggleBadChannel()");
    cw->unZoomButton->Connect("Clicked()", "GuiController", this, "UnZoom()");

    vw->can->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
        "GuiController",
        this,
        "ProcessCanvasEvent(Int_t,Int_t,Int_t,TObject*)"
    );
}

void GuiController::ToggleBadChannel()
{
    data->doDrawBadCh = cw->badChanelButton->IsDown();
    data->DrawBadCh();
}

void GuiController::SetCurrentCluster(int newCluster)
{
    cout << "new cluster: " << newCluster << endl;
    currentCluster = newCluster;
    data->currentCluster = currentCluster;
}

void GuiController::ClusterChanged(int i)
{
    int newCluster = cw->clusterEntry->GetNumber();
    if (newCluster == currentCluster) return;
    SetCurrentCluster(newCluster);

    data->DrawNewCluster();
}


void GuiController::UnZoom()
{
    data->ZoomProj(0, -1);
}

void GuiController::ZoomChanged()
{
    int zoomBin = 10*cw->zoomEntry->GetNumber();
    data->ZoomProj(currentPointIndex, zoomBin);
}

void GuiController::ProcessCanvasEvent(Int_t ev, Int_t x, Int_t y, TObject *selected)
{
    if (ev == 11) { // clicked
        if (!(selected->IsA() == TGraph::Class()
            || selected->IsA() == TH2F::Class()
            // || selected->IsA() == TLine::Class()
        )) return;
        TVirtualPad* pad = vw->can->GetClickSelectedPad();
        int padNo = pad->GetNumber();
        double xx = pad->AbsPixeltoX(x);
        double yy = pad->AbsPixeltoY(y);
        cout << "pad " << padNo << ": (" << xx << ", " << yy << ")" << endl;
        if (selected->IsA() == TGraph::Class() && padNo <=3) { // first row tgraph clicked
            TGraph *g = (TGraph*)gROOT->FindObject("g_dqdx");
            currentPointIndex = TMath::BinarySearch(g->GetN(), g->GetX(), xx);
            int zoomBin = 10*cw->zoomEntry->GetNumber();
            data->ZoomProj(currentPointIndex, zoomBin);
            data->DrawPoint(currentPointIndex);
        }
        else if (selected->IsA() == TH2F::Class()) { // 2nd row th2f clicked

        }
    }

}

void GuiController::HandleMenu(int id)
{
    // const char *filetypes[] = {"ROOT files", "*.root", 0, 0};
    switch (id) {
        case M_FILE_EXIT:
            gApplication->Terminate(0);
            break;
    }
}

TString GuiController::OpenDialog()
{
    const char *filetypes[] = {"ROOT files", "*.root", 0, 0};
    TString currentDir(gSystem->WorkingDirectory());
    static TString dir(".");
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    fi.fIniDir    = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), mw, kFDOpen, &fi);
    dir = fi.fIniDir;
    gSystem->cd(currentDir.Data());

    if (fi.fFilename) {
        // UnZoom();
        cout << "open file: " << fi.fFilename << endl;
        return fi.fFilename;
    }
    else {
        gApplication->Terminate(0);
    }
    return "";

}
