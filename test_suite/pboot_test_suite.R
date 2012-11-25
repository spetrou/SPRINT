library("RUnit")
library("boot")

data(city, package='boot')
data(gravity, package='boot')
data(nuclear, package='boot')
data(aircondit, package='boot')
data(discoveries)
data(trees)
for (nm in list.files("../inst/unitTests/pboot/", pattern = "\\.[Rr]$")){
  source(file.path("../inst/unitTests/pboot/", nm))
}

test.suite <- defineTestSuite("pboot", dirs = file.path("../inst/unitTests/pboot/"),testFileRegexp = '*.R')
library("sprint")
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
pterminate()
quit()
