# Deploys to shinyapps.io; see scripts/deploy.R for Posit Connect Cloud.
Sys.setenv(DEPLOY_SERVER = "shinyapps.io")
source("scripts/deploy.R")
