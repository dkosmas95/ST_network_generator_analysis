param numNodes;
param numNodesOrig;
param nPlans default 0;

set SOURCE;
set SINK;

set headNODES := 1..numNodes union SINK;				# non-(super)source nodes
set tailNODES := SOURCE union 1..numNodes;				# non-(super)sink nodes
param source symbolic in tailNODES;           			# (super)source
param sink symbolic in headNODES; 			   			# (super)sink

set trafficker within 1..numNodes;
set bottoms within 1..numNodes;
set victims within 1..numNodes;
set newTrafficker within 1..numNodes;
set replaceTrafficker within trafficker cross newTrafficker;
set newBottom within victims;
set replaceBottom within bottoms cross victims;

set ARCS within (tailNODES cross headNODES);      		# arcs
set SPLIT within (headNODES cross tailNODES);	    	# split nodes
set ARCSnewOut within (tailNODES cross headNODES);			# restructurable arcs
set ARCSnewIn within (tailNODES cross headNODES);

param budgetA;
param budgetD;								  
param cap {ARCS} >= 0;                        
param capnew {ARCSnewOut} >= 0;
param nodecap {SPLIT} >= 0;                   # capacities
param nodecapInc {newBottom} >= 0;
param cOut {ARCSnewOut} >= 0;
param cIn {ARCSnewIn} >= 0;

var hpi {0..nPlans, headNODES}; 							  # side of cut
var tpi {0..nPlans, tailNODES};
var theta {0..nPlans, ARCS union ARCSnewOut} >= 0;               # cut indicator
var nodetheta {0..nPlans, SPLIT} >= 0;

var flow {(i,j) in ARCS} >= 0, <= cap[i,j];      	# flows of existing arcs
var nodeflow {(i,j) in SPLIT} >=0;	# flows through nodes
var flownew {(i,j) in ARCSnewOut} >= 0;			 	# flows of restructured arcs

var y {SPLIT} binary;	  							# interdicted indicator
var zOut {(i,j) in ARCSnewOut} binary;	 				# restructured arcs Out
var zIn {(i,j) in ARCSnewIn} binary;				 	# restructured arcs In
var p {i in newBottom} binary;
var adjCost {i in trafficker} integer;
var intAdjCost {i in trafficker} >= 0;

param ybar{SPLIT} default 0;
param yBarM{0..nPlans, SPLIT} default 0;
param etaBar{0..nPlans} default 0;
param zbarOut {0..nPlans, ARCSnewOut} default 0;			# restructured indicators
param zbarIn {0..nPlans, ARCSnewIn} default 0;	
param pbar {0..nPlans, newBottom} default 0;

var deltaOut {0..nPlans, tailNODES} binary;
var deltaIn {0..nPlans, headNODES} binary;
var eta;

param r {tailNODES diff SOURCE};

param delayY {SPLIT};
param delayZOut {ARCSnewOut};
param delayZIn {ARCSnewIn};

var wOut {m in 0..nPlans, (i,j) in ARCSnewOut} binary; 							# indicator of ability to appear
var wIn {m in 0..nPlans, (i,j) in ARCSnewIn} binary;

param yopt {SPLIT} default 0;
param zoptOut {ARCSnewOut} default 0;
param zoptIn {ARCSnewIn} default 0;

## Master Problem
minimize objVal: eta;

subject to cutval {m in 0..nPlans}: 
  eta >= (sum {(i,j) in ARCS} cap[i,j]*theta[m,i,j] 
  + sum {(i,j) in ARCSnewOut} capnew[i,j]*theta[m,i,j] 
  + sum {(i,j) in SPLIT: i not in newBottom} nodecap[i,j]*nodetheta[m,i,j] 
  + sum {(i,j) in SPLIT: i in newBottom} (nodecap[i,j] + nodecapInc[i]*zbarOut[m,source,i]*wOut[m,source,i])*nodetheta[m,i,j]);

subject to sourceSide {m in 0..nPlans, (source,j) in ARCS}:
  hpi[m,j] + theta[m,source,j] >= 1;

subject to sourceSideOutNew {m in 0..nPlans, (source,j) in ARCSnewOut}:
  hpi[m,j] + theta[m,source,j] >= zbarOut[m,source,j] - 1 + wOut[m,source,j];

subject to sinkSide {m in 0..nPlans, (j, sink) in ARCS}:
  -1*tpi[m,j] + theta[m,j,sink] >= 0;
  
subject to arcInCut {m in 0..nPlans, (i,j) in ARCS: (i <> source and j <> sink)}:
  hpi[m,j] - tpi[m,i] + theta[m,i,j] >= 0;
  
subject to arcOutNewInCut {m in 0..nPlans, (i,j) in ARCSnewOut: (i <> source and j <> sink)}:
  hpi[m,j] - tpi[m,i] + theta[m,i,j] >= zbarOut[m,i,j] - 2 + wOut[m,i,j];

subject to arcInNewInCut {m in 0..nPlans, (i,j) in ARCSnewIn: (i <> source and j <> sink)}:
  hpi[m,j] - tpi[m,i] + theta[m,i,j] >= zbarIn[m,i,j] - 2 + wIn[m,i,j];
  
subject to nodeInCut {m in 0..nPlans, (i,j) in SPLIT}:
  tpi[m,j] - hpi[m,i] + nodetheta[m,i,j] >= -y[i,j];
  
subject to noIntNew {i in (numNodesOrig+1)..numNodes}:
  y[i,i] = 0;  
  
subject to knapsackA:
  sum {i in 1..numNodesOrig: i not in trafficker} r[i]*y[i,i] + sum{i in trafficker} intAdjCost[i] <= budgetA;

subject to intTrafCost1 {i in trafficker}:
  adjCost[i] >= r[i] - 3*sum{j in bottoms: (i,j) in ARCS}y[j,j] - sum{j in victims: (i,j) in ARCS}y[j,j];

subject to intTrafCost2 {i in trafficker}:
  adjCost[i] >= 4;
  
subject to trafIntMc1 {i in trafficker}:
  intAdjCost[i] <= r[i]*y[i,i];
  
subject to trafIntMc2 {i in trafficker}:
  intAdjCost[i] <= adjCost[i];

subject to trafIntMc3 {i in trafficker}:
  intAdjCost[i] >= r[i]*y[i,i] + adjCost[i] - r[i];

subject to reconstitutionOut {m in 1..nPlans, i in trafficker: sum{(i,h) in ARCSnewOut} zbarOut[m,i,h]>=1}:
  card({(i,h) in ARCS: h <> sink})*deltaOut[m,i] + sum{(i,j) in ARCSnewOut: zbarOut[m,i,j]=1} wOut[m,i,j] >= sum {(i,h) in ARCS: h <> sink} y[h,h];
   
subject to reconstitutionIn {m in 1..nPlans, j in victims: sum{(h,j) in ARCSnewIn} zbarIn[m,h,j]>=1}:
  card({(h,j) in ARCS: h <> source})*deltaIn[m,j] + sum{(i,j) in ARCSnewIn: i in trafficker and zbarIn[m,i,j]=1} wIn[m,i,j] >= sum {h in trafficker union bottoms: (h,j) in ARCS} y[h,h];

subject to fullOut {m in 1..nPlans, i in trafficker: sum{(i,h) in ARCSnewOut} zbarOut[m,i,h]>=1}:
  deltaOut[m,i] <= (sum{(i,j) in ARCSnewOut: zbarOut[m,i,j] = 1} wOut[m,i,j]) / card({(i,j) in ARCSnewOut: zbarOut[m,i,j]=1});

subject to fullIn {m in 1..nPlans, j in victims: sum{(h,j) in ARCSnewIn} zbarIn[m,h,j]>=1}:
  deltaIn[m,j] <= (sum{(i,j) in ARCSnewIn: i in trafficker and zbarIn[m,i,j] = 1} wIn[m,i,j]) / card({(i,j) in ARCSnewIn: zbarIn[m,i,j]=1});

subject to promotion {m in 1..nPlans, i in bottoms: sum{j in victims: (i,j) in replaceBottom} zbarOut[m,source,j] > 0}:
  sum {j in victims: (i,j) in replaceBottom and zbarOut[m,source,j] = 1} wOut[m,source,j] >= y[i,i];   
  
subject to newTraf {m in 1..nPlans, i in trafficker: sum{(i,l) in replaceTrafficker} zbarOut[m,source,l] > 0}:
  sum{j in newTrafficker: (i,j) in replaceTrafficker and zbarOut[m,source,j]=1} wOut[m,source,j] >= y[i,i];  
  
subject to passBotVicW {m in 1..nPlans, (i,j) in ARCSnewOut: i in bottoms}:
  wOut[m,i,j] >= zbarOut[m,i,j];
  
subject to passNewBotVicW {m in 1..nPlans, (i,j) in ARCSnewOut: i in newBottom}:
  wOut[m,i,j] >= zbarOut[m,i,j];

  
## Sub Problem
maximize outflow: 
  sum {(source,j) in ARCS} flow[source,j] 
  + sum {(source,j) in ARCSnewOut} flownew[source,j];

subject to inBalance {k in tailNODES diff {source}}:
  sum {(i,k) in ARCS} flow[i,k] + sum {(i,k) in ARCSnewOut} flownew[i,k] = nodeflow[k,k];

subject to outBalance {k in headNODES diff {sink}}:
  nodeflow[k,k] = sum {(k,j) in ARCS} flow[k,j] + sum {(k,j) in ARCSnewOut} flownew[k,j];
  
subject to nodeflowCap {(i,j) in SPLIT: i not in newBottom}:
  nodeflow[i,j] - nodecap[i,j]*(1 - ybar[i,j]) <= 0;
  
subject to nodeflowCapInc {(i,j) in SPLIT: i in newBottom}:
  nodeflow[i,j] - nodecap[i,j]*(1 - ybar[i,j]) - nodecapInc[i]*zOut[source,i] <= 0;
  
subject to flowCapNew {(i,j) in (ARCSnewOut inter ARCSnewIn)}:
  flownew[i,j] - capnew[i,j]*(zOut[i,j]+zIn[i,j]) <= 0;
  
subject to flowCapNewOut {(i,j) in ARCSnewOut diff ARCSnewIn}:
  flownew[i,j] - capnew[i,j]*zOut[i,j] <= 0;
  
subject to intNeighborBelow {i in trafficker: (i,sink) not in ARCS}:
  sum {(i,j) in ARCSnewOut} zOut[i,j] <= sum {(i,h) in ARCS: h <> sink} ybar[h,h];

subject to intNeighborAbove {i in victims: (source,i) not in ARCS}:
  sum {(j,i) in ARCSnewIn} zIn[j,i] <= sum {h in trafficker: (h,i) in ARCS} ybar[h,h];

subject to replaceIntTraf {(i,j) in replaceTrafficker}:
  zOut[source,j] <= ybar[i,i];

subject to replaceIntBottom {i in bottoms}:
  sum {j in victims: (i,j) in replaceBottom} zOut[source,j] <= ybar[i,i];

subject to noPromoteInt {i in newBottom}:
  zOut[source,i] <= 1 - ybar[i,i];

subject to newArcsT {i in trafficker}: 
  sum{(i,j) in ARCSnewOut} cOut[i,j]*zOut[i,j] + sum{(i,j) in ARCSnewIn} cIn[i,j]*zIn[i,j] + sum{k in newBottom: (i,k) in ARCS}cOut[source,k]*zOut[source,k] + sum{(k,j) in ARCSnewOut:(i,k) in ARCS} cOut[k,j]*zOut[k,j] + sum{(i,j) in replaceTrafficker}cOut[source,j]*zOut[source,j]  <= 8;

subject to newArcsTB {i in trafficker}:
  sum{j in newBottom: (i,j) in ARCS} zOut[source,j] <= 1;

subject to newArcsNewB {i in newBottom}:
  sum{j in victims: (i,j) in ARCSnewOut} zOut[i,j] <= card({j in victims: (i,j) in ARCSnewOut})*zOut[source,i]; 

subject to newArcsV {j in victims}:
  sum{i in 1..numNodes: (i,j) in ARCSnewIn} zIn[i,j] + sum{i in 1..numNodes: (i,j) in ARCSnewOut} zOut[i,j] <= 1;
   
subject to noRepeat {(i,j) in (ARCSnewOut inter ARCSnewIn)}:
  zIn[i,j]+zOut[i,j] <= 1;  

  
  