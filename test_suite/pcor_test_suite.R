library("RUnit")
library("ff")

for (nm in list.files("../inst/unitTests/pcor/", pattern = "\\.[Rr]$")){
  source(file.path("../inst/unitTests/pcor/", nm))
}
library("sprint")
test.suite <- defineTestSuite("pcor", dirs = file.path("../inst/unitTests/pcor/"),testFileRegexp = '*.R')

# === Set up finished ===

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)

pterminate()
quit()

