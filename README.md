
<b> <i> Mission Statement </b> </i>

Object-oriented C++/ROOT framework for making jet energy corrections and resolutions within the CMS Heavy Ions community. This repository measures dijet asymmetries and extrapolates out data to simulation differences using a pT-balance method. Currently usable for resolving jet energy and finding L2Residual corrections.
<br>

<i> Methodology and text-file convention follow CMS JERC tools (vendored in `external/jetmet*`)
<br>
see [CMS JERC](https://cms-jerc.web.cern.ch/) for official recommendations and definitions </i>

<h3> Workflow </h3>

||||
| :-: | :-: | --- |
| Step 1 | `runAsymmetry` (`make_condor.sh`) | find dijet asymmetries (with HTCondor) |
| <i> Step 1.5 </i> | `batch_hadd.sh` | <i> combine HTCondor output </i> |
| Step 2 | `runCalibration` | extract L2Residuals or JER SF from dijet asymmetries |
| <i> Step 2.5 </i> | `runResponse` | <i> extract JES and JER from MC dijet asymmetries </i> |
| Step 3 | `runTextFile` | make plain text files from calibration |
| <i> n/a </i> | `runPlotting` | make png files from any output |

<h2> Usage </h2>

Copy and paste code blocks for generic workflow
- -config is required everywhere, cfg/default.toml is a template config file
    - see cfg/2024ppref.toml for working example
- consider | as OR
<br><br>

<strong> Build </strong>

<i> SCRAM </i>

```bash
cd CMSSW_X_Y_Z/src
cmsenv
mkdir -p Analysis
git clone git@github.com:thenicholasbarnett/HIN-L2Residuals-JER.git Analysis/HIN-L2Residuals-JER
cd Analysis/HIN-L2Residuals-JER
scram b
```

<i> CMake </i>

```bash
git clone git@github.com:thenicholasbarnett/HIN-L2Residuals-JER.git
cd HIN-L2Residuals-JER
cmake -B build
cmake --build build
```

<i> Rebuild </i>

```bash
cmsenv && scram b clean && scram b
```

```bash
rm -rf build && cmake -B build && cmake --build build
```

<strong> Step 1 </strong>

<i> HTCondor </i>

```bash 
bash condor/make_condor.sh \
  -output /path/to/output/dir \
  -filelists path/to/filelist_triggered.txt \
             path/to/filelist_non-triggered.txt \
             path/to/filelist_MC.txt \
  -tag subset \
  -config cfg/configuration.toml
```

```bash 
bash condor/make_condor.sh \
  -output /path/to/output/dir \
  -filelists path/to/filelist/dir \
  -config cfg/configuration.toml
```

```bash
bash condor/batch_hadd.sh \
  -output /path/to/output/merged-output.root \
  -input "/path/to/output/dir/*.root" \
  -batchsize 10 \
  -njobs 2
```

<i> Local </i>

```bash
./build/bin/runAsymmetry \
  -input path/to/HiForest_mc|triggered|non-triggered.root \
  -output path/to/output_mc|triggered|non-triggered.root \
  -mode mc|triggered|non-triggered \
  -config cfg/configuration.toml
```

<i> JER SF closure (MC only) </i>

```bash
./build/bin/runAsymmetry \
  -input path/to/HiForest_mc.root \
  -output path/to/residuals-closure_mc.root \
  -mode mc \
  -calibration jer \
  -closure true \
  -config cfg/configuration.toml
```

`-calibration jer -closure true` JER-smears MC jets (hybrid method,
`JetSmearing.h`) before histogramming, using
`jer_closure.resolution_files`/`scale_factor_files` from the config, to
produce the smeared-MC half of the JER SF closure check fed into
`runCalibration -mode jer` (see Step 2 below). Only valid with `-mode mc`;
any other `-calibration`/`-closure` combination is a no-op.

<strong> Step 2 </strong>

<i> Residuals </i>

```bash
./build/bin/runCalibration \
  -data path/to/asymmetry-triggered.root \
  -mc path/to/asymmetry-mc.root \
  -output data/root/second/residuals-triggered.root \
  -mode jec \
  -config cfg/residuals-configuration.toml
```

<i> Response </i>

```bash
./build/bin/runResponse \
  -input path/to/asymmetry-MC.root \
  -output path/to/response-MC.root \
  -config cfg/residuals_configuration.toml
```

<i> Resolution </i>

```bash
./build/bin/runCalibration \
  -data path/to/residuals-closure_triggered.root \
  -mc path/to/residuals-closure_mc.root \
  -output data/root/second/resolution_triggered.root \
  -mode jer \
  -config cfg/resolutions-configuration.toml
```

<strong> Step 3 </strong>

<i> Single Dataset Type </i>

```bash
./build/bin/runTextFile \
  -triggered|non-triggered path/to/residuals|resolutions_triggered|non-triggered.root \
  -output path/to/corrections|resolutions.root \
  -mode jec|jer \
  -config cfg/configuration.toml
```

<i> Multiple Dataset Types </i>

```bash
./build/bin/runTextFile \
  -triggered path/to/residuals|resolutions_triggered.root \
  -nontriggered path/to/residuals|resolutions_non-triggered.root \
  -output path/to/corrections|resolutions.root \
  -mode jec|jer \
  -config cfg/configuration.toml
```

<strong> Plotting </strong>

```bash
./build/bin/runPlotting \
  -input path/to/input.root \
  -outdir path/to/output/dir \
#  -flags "all" \
  -config cfg/configuration.toml \
```

<strong> Configuration File </strong>

```bash
./build/bin/runPlotting args.config
```

```
input = path/to/input.root
outdir = path/to/output/dir
config = cfg/configuration.toml
```


<h2> Command Line Keys </h2>

Execute binaries by providing CL keys with values (-key = value)
<br>
Vendored `CommandLine.h` from JetMET, can be found in external/jetmet

<table>
  <thead>
    <tr>
      <th align="center">Key</th>
      <th align="center">Used by</th>
      <th align="center">Values</th>
      <th>Meaning</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center"><code>config</code></td>
      <td align="center">all compiled binaries<br><code>make_condor.sh</code></td>
      <td align="center">TOML path</td>
      <td>Analysis configuration. Always required.</td>
    </tr>
    <tr>
      <td align="center"><code>input</code></td>
      <td align="center"><code>runAsymmetry</code><br><code>runResponse</code><br><code>runPlotting</code></td>
      <td align="center">ROOT path</td>
      <td>Input ROOT file.</td>
    </tr>
    <tr>
      <td align="center"><code>output</code></td>
      <td align="center"><code>runAsymmetry</code><br><code>runCalibration</code><br><code>runResponse</code><br><code>runTextFile</code><br><code>make_condor.sh</code></td>
      <td align="center">path</td>
      <td>Output ROOT file or output directory.</td>
    </tr>
    <tr>
      <td align="center"><code>mode</code></td>
      <td align="center"><code>runAsymmetry</code></td>
      <td align="center"><code>triggered</code><br><code>non-triggered</code><br><code>mc</code></td>
      <td>Selects data/MC event handling.</td>
    </tr>
    <tr>
      <td align="center"><code>mode</code></td>
      <td align="center"><code>runCalibration</code><br><code>runTextFile</code></td>
      <td align="center"><code>jec</code><br><code>jer</code></td>
      <td>Selects mean-derived JEC or width-derived JER output. Required for <code>runTextFile</code> whenever <code>triggered</code>/<code>nontriggered</code> is used.</td>
    </tr>
    <tr>
      <td align="center"><code>maxevents</code></td>
      <td align="center"><code>runAsymmetry</code></td>
      <td align="center">integer</td>
      <td>Limits processed events. Default: all.</td>
    </tr>
    <tr>
      <td align="center"><code>data</code></td>
      <td align="center"><code>runCalibration</code></td>
      <td align="center">ROOT path</td>
      <td>Hadded Step 1 data file.</td>
    </tr>
    <tr>
      <td align="center"><code>mc</code></td>
      <td align="center"><code>runCalibration</code></td>
      <td align="center">ROOT path</td>
      <td>Hadded Step 1 MC file.</td>
    </tr>
    <tr>
      <td align="center"><code>triggered</code></td>
      <td align="center"><code>runTextFile</code></td>
      <td align="center">ROOT path</td>
      <td>Triggered Step 2 residual file.</td>
    </tr>
    <tr>
      <td align="center"><code>nontriggered</code></td>
      <td align="center"><code>runTextFile</code></td>
      <td align="center">ROOT path</td>
      <td>Non-triggered Step 2 residual file.</td>
    </tr>
    <tr>
      <td align="center"><code>resolution</code></td>
      <td align="center"><code>runTextFile</code></td>
      <td align="center">ROOT path</td>
      <td><code>runResponse</code> output used for JER pT-resolution text.</td>
    </tr>
    <tr>
      <td align="center"><code>tag</code></td>
      <td align="center"><code>runTextFile</code><br><code>runPlotting</code><br><code>make_condor.sh</code></td>
      <td align="center">plain name</td>
      <td>Output label. No slashes.</td>
    </tr>
    <tr>
      <td align="center"><code>method</code></td>
      <td align="center"><code>runTextFile</code></td>
      <td align="center"><code>gauss</code><br><code>doubleGauss</code><br><code>trunc90</code><br><code>trunc95</code></td>
      <td>Step 3 fit source. Default comes from TOML.</td>
    </tr>
    <tr>
      <td align="center"><code>norm</code></td>
      <td align="center"><code>runTextFile</code></td>
      <td align="center"><code>true</code><br><code>false</code></td>
      <td>Uses normalized intercepts when true. Default: true.</td>
    </tr>
    <tr>
      <td align="center"><code>outdir</code></td>
      <td align="center"><code>runPlotting</code></td>
      <td align="center">directory</td>
      <td>Plot output directory.</td>
    </tr>
    <tr>
      <td align="center"><code>flags</code></td>
      <td align="center"><code>runPlotting</code></td>
      <td align="center"><code>all</code><br>or keywords</td>
      <td>Plot families to draw.</td>
    </tr>
    <tr>
      <td align="center"><code>closure</code></td>
      <td align="center"><code>runPlotting</code><br><code>runAsymmetry</code></td>
      <td align="center"><code>true</code><br><code>false</code></td>
      <td><code>runPlotting</code>: fixes final-correction y range for closure checks.<br><code>runAsymmetry</code>: combined with <code>calibration jer</code>, JER-smears MC jets before histogramming. Default: false.</td>
    </tr>
    <tr>
      <td align="center"><code>calibration</code></td>
      <td align="center"><code>runPlotting</code><br><code>runAsymmetry</code></td>
      <td align="center"><code>JEC</code>/<code>jec</code><br><code>JER</code>/<code>jer</code></td>
      <td><code>runPlotting</code>: selects plotted calibration objects.<br><code>runAsymmetry</code>: only meaningful combined with <code>-mode mc -closure true</code> (see above). Default: JEC/jec.</td>
    </tr>
    <tr>
      <td align="center"><code>sample</code></td>
      <td align="center"><code>runPlotting</code></td>
      <td align="center"><code>true</code><br><code>false</code></td>
      <td>Stops after 1000 plots or 60s, whichever comes first. Default: false.</td>
    </tr>
    <tr>
      <td align="center"><code>alltxt</code></td>
      <td align="center"><code>make_condor.sh</code></td>
      <td align="center">flag</td>
      <td>Uses every filelist in <code>data/txt/</code>.</td>
    </tr>
    <tr>
      <td align="center"><code>filelists</code></td>
      <td align="center"><code>make_condor.sh</code></td>
      <td align="center">paths</td>
      <td>Uses selected filelists or directories.</td>
    </tr>
    <tr>
      <td align="center"><code>nosubmit</code></td>
      <td align="center"><code>make_condor.sh</code></td>
      <td align="center">flag</td>
      <td>Writes Condor files without submitting.</td>
    </tr>
  </tbody>
</table>

Plot flags:

<p align="center"><code>etasym</code> · <code>methods</code> · <code>finals</code> · <code>normcomp</code> · <code>adist</code> · <code>roverlay</code> · <code>alpha</code> · <code>ptfit</code> · <code>kinematics</code> · <code>event</code> · <code>response</code></p>
