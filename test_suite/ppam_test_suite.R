library("RUnit")
library("ff")
library("cluster")

for (nm in list.files("../inst/unitTests/ppam/", pattern = "\\.[Rr]$")){
  source(file.path("../inst/unitTests/ppam/", nm))
}
library("sprint")
test.suite <- defineTestSuite("ppam", dirs = file.path("../inst/unitTests/ppam/"),testFileRegexp = '*.R')

# === Set up finished ===

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)

pterminate()
quit()

