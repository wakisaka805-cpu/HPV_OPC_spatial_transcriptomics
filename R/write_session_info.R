source(file.path("R", "config.R"))
capture.output(sessionInfo(), file = file.path(REPO_ROOT, "sessionInfo.txt"))
message("Wrote sessionInfo.txt")
