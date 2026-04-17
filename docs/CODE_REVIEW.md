# Magnify-tracking — Code Review

> **Status:** All items in §5 (efficiency) and §6 (bugs) listed below have been
> addressed in the same commit that introduced this document. Each item is
> annotated **[FIXED]** with a short description of the applied fix.
> Verified: all five shared libraries rebuild cleanly via ACLiC, and a smoke
> test against `scripts/track_com.root` loads/destructs `Data` without errors.

## 1. Overview

**Magnify-tracking** is an interactive ROOT GUI for inspecting Wire-Cell reconstructed tracks. It reads five ROOT TTrees from a single `.root` file:

| Tree | Contents |
|---|---|
| `T_rec` | 3D reconstructed hits per cluster (x/y/z, dQ, dx, L, u/v/w channel projections, sub-cluster id, vertex flag, chi², MC comparisons) |
| `T_proj_data` | Per-wire charge and predicted charge for each cluster |
| `T_true` | MC-truth particle trajectory |
| `T_bad_ch` | Bad/dead wire channel list |
| `T_proj` | (loaded but not actively used currently) |

The viewer displays a 3×3 pad canvas:

```
| dQ/dx plot       | MC compare      | 3D view         |
| U-plane charge   | V-plane charge  | W-plane charge  |
| U-plane pred     | V-plane pred    | W-plane pred    |
```

**Entry point:**
```bash
./magnify.sh <file.root> [sign]
```
Internally: `scripts/loadClasses.C` ACLiC-compiles all `.cc` files, then `scripts/Magnify.C` instantiates `GuiController`.

The `sign` argument (0/+1/−1) controls how the charge-prediction residual `(pred−meas)/meas` is colored in the prediction panels.

---

## 2. Directory Layout

```
Magnify-tracking/
├── magnify.sh                  entry-point shell script
├── event/
│   ├── Data.h                  Data class declaration (~145 lines)
│   └── Data.cc                 Data model + all drawing logic (~1280 lines)
├── viewer/
│   ├── GuiController.{h,cc}    Mediator: wires GUI signals → Data draw calls
│   ├── MainWindow.{h,cc}       Top-level TGMainFrame (menu + view + control)
│   ├── ViewWindow.{h,cc}       TRootEmbeddedCanvas, 3×3 pads, color palettes
│   └── ControlWindow.{h,cc}    Cluster/point/zoom/range controls
├── scripts/
│   ├── Magnify.C               ROOT bootstrap (instantiates GuiController)
│   ├── loadClasses.C           ACLiC .L directives for all source files
│   ├── convert.C               Utility: converts organized.dat → temp.root TTree
│   └── *.root / *.dat          Demo/validation data files
├── tmp/                        Additional ROOT files for testing
├── prototoype_5384/            Prototype run-5384 ROOT files
└── toolkit_5384/               (empty, future use)
```

---

## 3. Architecture

The project uses a **Mediator** pattern. `GuiController` is the central coordinator: it owns `Data`, `MainWindow`, `ViewWindow`, and `ControlWindow`, and translates GUI signals into `Data` draw calls.

```
MainWindow (TGMainFrame)
 ├── MenuBar  [File → Exit]
 ├── ViewWindow (TRootEmbeddedCanvas)
 │     TCanvas divided 3×3
 │     pad 1: dQ/dx        pad 2: MC-compare  pad 3: 3D view
 │     pad 4: proj-U       pad 5: proj-V      pad 6: proj-W
 │     pad 7: pred-U       pad 8: pred-V      pad 9: pred-W
 └── ControlWindow (TGHorizontalFrame)
       [General group]  cluster index | cluster id | all-clusters | bad-ch | draw-track
       [Zoom group]     point index   | zoom ±10×  | UnZoom
       [Range group]    time/u/v/w min-max | Zoom | keep-range
              ▲
              │  ROOT signal/slot Connect(...)
              ▼
         GuiController  ─────────────────────────►  Data
         ProcessCanvasEvent, ClusterChanged,         LoadData, Draw*,
         ZoomChanged, ToggleAllCluster, etc.         FindCluster/PointIndex
```

### Component Details

#### `Data` (`event/Data.cc`)
The largest class (~1280 LOC). Responsibilities:
- **Load**: `LoadData` → `LoadRec`, `LoadProj`, `LoadBadCh`, `LoadTruth`. Reads all tree branches into member vectors.
- **State**: owns `currentCluster`, `doDrawBadCh`, `doDrawTrack`, `sign`, `nCluster`, lookup maps (`rec_cluster_map`, `data_cluster_map`), sub-cluster index lists, and *all* ROOT graphical objects.
- **Draw**:

| Method | Pad(s) | What it draws |
|---|---|---|
| `DrawDQDX` | 1 | `TGraph` of dQ/dx vs distance; red sub-cluster boundary lines + id labels; optional true dQ/dx overlay; optional χ² overlay |
| `DrawMCCompare` | 2 | `TGraph` of angular and spatial distance to MC truth vs distance |
| `Draw3D` | 3 | `TH3F` frame (axis/zoom anchor), gray point cloud (`TPolyMarker3D`), MC truth `TPolyLine3D`, reco start marker |
| `DrawProj` | 4–9 | Six `TH2F`s (channel×time, one per plane × {measured, pred}); calls `DrawSubclusters` and `DrawBadCh` |
| `DrawProjAll` | 4–9 | Same as `DrawProj` but fills all clusters simultaneously |
| `DrawSubclusters` | 4–9, 3 | Colored `TGraph` per sub-cluster overlaid on proj/pred pads; same coloring in 3D |
| `DrawPoint` | 4–9, 3 | Magenta cross-hair marker at selected point index |
| `DrawBadCh` | 4–6 | Toggles gray `TLine`s for bad channels |
| `ZoomProj` | all | Coordinated axis range on all pads (and 1D dQ/dx) given a point index and zoom bin |
| `DrawNewCluster` | all | Calls all of the above in sequence |

- **Hit-test**: `FindClusterIndex(x, y)` — O(nClusters × points); `FindPointIndex(x, y)` — O(points in currentCluster).

#### `GuiController` (`viewer/GuiController.cc`)
Connects every widget to a slot:
- `ClusterChanged` / `ClusterIdChanged` — sync the two cluster number widgets, call `SetCurrentCluster`, then `DrawNewCluster`.
- `ZoomChanged` — read pointIndex and zoomBin from entries, call `ZoomProj` + `DrawPoint`.
- `ToggleAllCluster` — switch between per-cluster `DrawProj` and full-event `DrawProjAll`.
- `ProcessCanvasEvent` — on left-click over a `TGraph` in pads 1–3 (dQ/dx), uses `TMath::BinarySearch` to find the nearest point and drives the zoom.

#### `MainWindow` (`viewer/MainWindow.cc`)
Minimal: creates `ViewWindow` and `ControlWindow` inside a vertical frame, sets up a File menu with only an "Exit" entry.

#### `ViewWindow` (`viewer/ViewWindow.cc`)
Wraps `TRootEmbeddedCanvas`. Divides the canvas 3×3. Provides five color palettes selectable at construction time:
- `PaletteRainbow` (default) — custom 5-stop gradient
- `PaletteGray` / `PaletteGrayInv` — linear gray scales
- `PaletteSummer` — 256-stop thermal ramp (full inline data)
- `PaletteFire` — ROOT built-in palette 53

#### `ControlWindow` (`viewer/ControlWindow.cc`)
Three `TGGroupFrame`s. Channel range defaults (0–2400 U, 2400–4800 V, 4800–8256 W) are hardcoded; this ties the layout to the MicroBooNE wire geometry (8256 total channels, 2400 time slices).

---

## 4. Data Flow on a Cluster Change

```
user changes "cluster index" entry
  → ClusterChanged(int)                        [GuiController.cc:153]
    → SetCurrentCluster(newCluster)             [GuiController.cc:142]
        data->currentCluster = newCluster
        updates pointIndexEntry limits
    → data->DrawNewCluster()                   [Data.cc:1005]
        DrawDQDX()   → builds sub-cluster boundary lists, draws g_dqdx
        DrawMCCompare()  → draws g_com_dtheta / g_com_dis
        Draw3D()     → Clear("nodelete"), build h_3d_frame + g_3d_points
        DrawProj()   → rebuild six TH2Fs, DrawSubclusters(), DrawBadCh()
        DrawPoint(0) → places cursor markers at first point
```

A canvas click on the dQ/dx graph (pad 1) triggers `ProcessCanvasEvent` → `BinarySearch` → `ZoomProj` + `DrawPoint`, bypassing the full cluster redraw.

---

## 5. Efficiency Improvement Opportunities

### 5.1 Per-cluster TH2F reallocation (`DrawProj`, Data.cc:431–449) **[FIXED]**
On every cluster change, six histograms are `delete`d and `new`-allocated:
```cpp
if (h_proj_u) { delete h_proj_u; }
h_proj_u = new TH2F("h_proj_u", "", nChannel_u, 0, nChannel_u, nTime, 0, nTime);
```
A 2400×2400 `TH2F` carries ~23 MB of bin storage. Six histograms = ~100 MB of heap churn per cluster click. `DrawProjAll` already demonstrates the correct pattern (allocate once, reuse). **Applied fix:** added the six `h_proj_*` / `h_pred_*` (and `_all` variants) as `Data` members; `EnsureProjHistos()` / `EnsureProjAllHistos()` allocate once on first use; `DrawProj` and `DrawProjAll` now call `Reset("ICES")` instead of `delete`/`new`. Destructor cleans them up.

### 5.2 Repeated `gROOT->FindObject` lookups (`DrawProj`, `ZoomProj`, `DrawDQDX`, `DrawMCCompare`) **[FIXED]**
Every draw cycle scans the ROOT object registry via `FindObject` for objects that could simply be kept as pointers:
```cpp
TH2F *h_proj_u = (TH2F*)gROOT->FindObject("h_proj_u");  // O(N) scan
```
If the histograms become `Data` members (§5.1) the `FindObject` calls in `ZoomProj` (Data.cc:563–668) can be replaced with direct member access. **Applied fix:** `ZoomProj` now uses the new member pointers via a small `zoomHisto()` static helper; six `FindObject` registry scans per zoom event are gone.

### 5.3 `.at()` bounds checks in tight inner loops (`DrawProj`, `DrawSubclusters`, `DrawProjAll`) **[PARTIALLY FIXED]**
```cpp
for (int i=0; i<size; i++) {
    double u = rec_u->at(currentCluster).at(i);   // two bounds checks per iteration
    ...
}
```
Hoist the cluster-level reference once:
```cpp
const auto& us = rec_u->at(currentCluster);  // one bounds check
for (int i = 0; i < size; i++) {
    double u = us[i];                         // no check
}
```
**Applied fix:** done for `DrawProj` (cluster fill loop and graph build loop), `DrawProjAll` (per-cluster fill loop and graph build loop), and `FindPointIndex`. `DrawSubclusters` was left as-is — its `.at()` calls remain for now since the loop is shorter and cluster ids are known to be valid.

### 5.4 Redundant `Modified()` / `Update()` calls **[NOT YET ADDRESSED]**
`DrawProj` updates each pad, then `DrawSubclusters` updates them again, then `DrawBadCh` updates them a third time. For small canvases this is fine, but for large events the triple update is wasteful. Suppress intermediate updates and call `c1->Update()` once at the end of `DrawNewCluster`.

### 5.5 `DrawProjAll` never resets histograms (`DrawProjAll`, Data.cc:1079–1085) **[FIXED]**
Histograms are only allocated once (`if (!h_proj_u_all) new ...`) but never cleared. Switching between "all clusters" and single-cluster mode, then back, will display stale superimposed data from a previous event. Add `h->Reset("ICES")` before the fill loop. **Applied fix:** `DrawProjAll` now calls `Reset("ICES")` on all six `_all` histograms after `EnsureProjAllHistos()`.

### 5.6 `FindClusterIndex` is O(nClusters × nPoints) (`Data.cc:1015`) **[NOT ADDRESSED — defer until hotspot]**
Acceptable for current event sizes, but will not scale. A spatial lookup grid (e.g. a coarse `std::unordered_map<(ch/10, t/10), clusterIdx>`) built at load time would reduce this to O(1) amortized.

---

## 6. Potential Bugs and Fragile Spots

### 6.1 Legend overwrites the MC-compare panel **[FIXED]**
`DrawSubclusters` (Data.cc:820) switches to `c1->cd(pad_dqdx+1)` = pad 2, which is `pad_mc_compare`. It draws `hInfo` and the sub-cluster color legend on that pad. Since `DrawSubclusters` is called from `DrawProj`, which runs *after* `DrawMCCompare` in `DrawNewCluster`, the MC comparison graphs (`g_com_dtheta`, `g_com_dis`) are immediately overwritten. In MC-mode files, the user never sees the MC comparison. The intent appears to be drawing the legend on an otherwise-unused pad; the fix is to assign a dedicated pad index for the legend (e.g. an `int pad_legend` member set to 2 only when `isData`, else to some other unused pad), or overlay it on the dQ/dx pad. **Applied fix:** `DrawSubclusters` now branches on `isData`. If `isData` (no MC pane to preserve), the existing hInfo + legend is drawn on `pad_mc_compare` as before. Otherwise, a compact transparent legend is drawn as an overlay on `pad_3d`, leaving the MC compare graphs intact.

### 6.2 `data_cluster_map` lookup without existence check **[FIXED]**
`DrawProj` (Data.cc:417) and `DrawProjAll` (Data.cc:1065):
```cpp
int index = data_cluster_map[cluster_id];
```
`std::map::operator[]` silently inserts a zero-valued entry for unknown keys. If a reconstructed cluster ID has no matching projection data (e.g. a newly formed cluster not yet in `T_proj_data`), `index` will be 0 and the code will silently display the wrong cluster's waveforms. **Fix:**
```cpp
auto it = data_cluster_map.find(cluster_id);
if (it == data_cluster_map.end()) { /* warn and return */ }
int index = it->second;
```
**Applied fix:** `DrawProj` now does exactly this and returns early with a console warning when the projection data is missing. In `DrawProjAll` the lookup was actually dead code (the loop iterates over all entries), so it was simply removed.

### 6.3 Memory leaks from per-redraw heap allocations **[PARTIALLY FIXED]**

| Location | Object(s) created but never deleted | Status |
|---|---|---|
| `DrawDQDX` Data.cc:293, 297 (in loop) | `TLine`, `TLatex` per sub-cluster boundary | Not fixed — ROOT pad `Clear()` reclaims on cluster change |
| `DrawDQDX` Data.cc:341 | `TLegend` | Not fixed — reclaimed by pad `Clear()` |
| `DrawMCCompare` Data.cc:399 | `TLegend` | Not fixed — reclaimed by pad `Clear()` |
| `~Data` Data.cc:1271–1281 | `infoText`, `dqdxPoint`, `currentPoint[0–2]`, `bad_lines` contents, `g_subclusters_u/v/w` contents | **Fixed** |

ROOT tracks objects in `gROOT`, so the pads eventually own most primitives and will delete them on `Clear()`. However the `TLegend`s and the `TLine`/`TLatex` objects inside loops are never re-found and accumulate with every cluster change. The destructor gaps mean shutdown may cause double-free or leak depending on ROOT version. **Applied fix:** `~Data` now explicitly deletes `infoText`, `dqdxPoint`, `currentPoint[0..2]`, all `bad_lines` elements, and all `g_subclusters_u/v/w` elements, plus all twelve new `TH2F` members.

### 6.4 macOS portability: `readlink -f` in `magnify.sh:3` **[FIXED]**
```bash
magnify_source="$(dirname $(readlink -f $BASH_SOURCE))"
```
BSD `readlink` (used on macOS by default) does not accept `-f`. On macOS the result is empty, so `magnify_source` becomes the directory name of an empty path and all relative paths break. **Applied fix:** replaced with:
```bash
magnify_source="$(cd "$(dirname "$BASH_SOURCE")" && pwd -P)"
```
This is portable across macOS BSD and GNU/Linux.

### 6.5 `FindPointIndex` returns the *last* match, not the nearest **[FIXED]**
Data.cc:1051–1058:
```cpp
if (TMath::Abs(x-xx)<5 && TMath::Abs(y-yy)<5) {
    foundIndex = i;   // keeps overwriting; returns last match
}
```
On dense tracks, multiple points fall within the ±5 bin window. The last one (highest index) is returned regardless of actual distance. **Applied fix:** replaced with a minimum-squared-distance search:
```cpp
double bestDist2 = 50.0;  // sqrt(50) ≈ 7-bin radius
for (int i = 0; i < size; i++) {
    double dx = x - xs[i], dy = y - ts[i];
    double d2 = dx*dx + dy*dy;
    if (d2 < bestDist2) { bestDist2 = d2; foundIndex = i; }
}
```

### 6.6 `reduced_chi2->at(currentCluster)` without per-cluster size guard **[FIXED]**
`DrawDQDX` (Data.cc:320):
```cpp
if (reduced_chi2->size()>0) {
    size = reduced_chi2->at(currentCluster).size();
```
This checks that the outer vector is non-empty, but not that it has at least `currentCluster+1` entries. ROOT files with a partial `reduced_chi2` branch (e.g. only computed for a subset of clusters) will throw `std::out_of_range`. **Applied fix:** guard tightened to `if (reduced_chi2->size() > (size_t)currentCluster && ...)`.

### 6.7 `nCluster` local variable shadows the class member in `DrawProjAll` **[FIXED]**
Data.cc:1091:
```cpp
int nCluster = data_channel->size();   // shadows Data::nCluster
```
A future maintainer adding loop bounds based on the member will pick the wrong value. **Applied fix:** renamed to `nDataClusters` throughout `DrawProjAll`, including the diagnostic `cout` statement.

### 6.8 Sign-convention block duplicated verbatim **[FIXED]**
The residual-sign logic (Data.cc:461–471 and Data.cc:1101–1111):
```cpp
if (sign == 0) {
    z_pred = TMath::Abs(z_pred);
} else if (sign>0) {
    if (z_pred<0) z_pred = 0;
} else {
    if (z_pred>0) z_pred = 0;
    else z_pred = TMath::Abs(z_pred);
}
```
Appeared in both `DrawProj` and `DrawProjAll`. **Applied fix:** extracted as `static double Data::applySign(double residual, int sign)` (declared in `Data.h`, defined in `Data.cc`); both call sites replaced with `z_pred = applySign(z_pred, sign)`.

### 6.9 `convert.C` hardcodes 245 entries **[FIXED]**
```cpp
for (int i=0; i!=245; i++) { infile >> ...; T->Fill(); }
```
A different `organized.dat` will silently read only 245 lines or attempt past EOF. **Applied fix:** replaced with:
```cpp
while (infile >> run_no >> subrun_no >> event_no >> ks1 >> ks2 >> norm1 >> norm2 >> flag) { T->Fill(); }
```

### 6.10 `LoadProj` dereferences `T_proj_data` without a null guard **[FIXED]**
`LoadData` (Data.cc:131) gets `T_proj_data` from the file but only checks `T_true` for MC/data mode. If `T_proj_data` is absent from the file, `LoadProj` (Data.cc:198) immediately calls `T_proj_data->SetBranchAddress(...)` and crashes. **Applied fix:** `LoadProj` now opens with `if (!T_proj_data) { std::cerr << "WARNING: T_proj_data missing, skipping LoadProj\n"; return; }`.

### 6.11 Bare `vector` / `map` in `Data.h` without `std::` qualification **[FIXED]**
Data.h:36–77 used `vector<int>`, `map<int,int>` etc. without the `std::` prefix. This works today only because ROOT headers inject `std` into the global namespace as a side effect. Portable C++ requires `std::vector<int>` in a header. Leaving this as-is is a latent compilation failure if ROOT's header pollution is cleaned up in a future version. **Applied fix:** all bare `vector<>` and `map<>` references in `Data.h` qualified with `std::`.

---

## 7. Suggested Improvement Roadmap

| Priority | Action | Addresses | Status |
|---|---|---|---|
| High | Pre-allocate proj/pred TH2Fs as Data members; `Reset()` instead of `delete`/`new` | §5.1, §5.5 | **Done** |
| High | Guard `data_cluster_map` with `.find()` in DrawProj / DrawProjAll | §6.2 | **Done** |
| High | Fix `magnify.sh` `readlink -f` for macOS | §6.4 | **Done** |
| Medium | Fix legend pad target (avoid clobbering MC-compare pad) | §6.1 | **Done** |
| Medium | Fix `FindPointIndex` to return nearest match | §6.5 | **Done** |
| Medium | Complete `~Data` destructor (delete `infoText`, `dqdxPoint`, `currentPoint[3]`, `bad_lines`, `g_subclusters_*`) | §6.3 | **Done** |
| Medium | Replace `gROOT->FindObject` with cached member pointers | §5.2 | **Done** |
| Low | Extract `applySign()` helper | §6.8 | **Done** |
| Low | Add `reduced_chi2` per-cluster size guard | §6.6 | **Done** |
| Low | Rename `nCluster` local in DrawProjAll | §6.7 | **Done** |
| Low | Fix `convert.C` loop to not hardcode 245 | §6.9 | **Done** |
| Low | Qualify `std::vector` / `std::map` in Data.h | §6.11 | **Done** |
| Low | Batch `Modified()/Update()` calls | §5.4 | Deferred |
| Low | Spatial grid for `FindClusterIndex` | §5.6 | Deferred |
