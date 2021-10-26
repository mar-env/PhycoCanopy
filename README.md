# PhycoCanopy
R Shiny code for PhycoCanopy macroalgal simulation

The libraries for shiny and REdaS should be loaded first in R or RStudio
library(shiny) 
library(REdaS) 

The following code will load the PhycoCanopy app when executed in R
runGitHub("PhycoCanopy","mar-env",ref="main")

The programme will not calculate if depth and tidal settings cause the water depth over the point where the algae are attached to exceed 10 m.

Be (relatively) patient if you choose to simulate periods longer than 1 day. It can take 10 seconds to recalculate the response to changed parameters for a 14 day simulation (3.1 GHz clock speed computer).

