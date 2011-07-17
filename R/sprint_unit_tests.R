##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

sprint.test_ppam <- function() {

   if( !require("RUnit", quietly=TRUE) ) { 
        warning("Unit tests not run for ppam - failed to load package RUnit")
        return()
    }

    if( !require("sprint", quietly=TRUE) ) {
        warning("Unit tests not run for ppam - failed to load package sprint")
        return()
    }

   if( !require("cluster", quietly=TRUE) ) {
     warning("Unit tests not run for ppam - failed to load package cluster")
     return()
   }

   if( !require("ff", quietly=TRUE) ) {
     warning("Unit tests not run for ppam - failed to load package ff")
     return()
   }

   # Find path for test scripts
   path <- system.file(package="sprint", "unitTests", "ppam")

   # Create test suite and execute it
   testSuite <- defineTestSuite(name="ppam", dirs=path)
   testData  <- runTestSuite(testSuite)
  
   # Print out the results
   printTextProtocol(testData, showDetails=TRUE)

}

sprint.test_pmaxT <- function() {

    if( !require("RUnit", quietly=TRUE) ) { 
        warning("Unit tests not run for pmaxT - failed to load package RUnit")
        return()
    }

    if( !require("sprint", quietly=TRUE) ) {
        warning("Uunit tests not run for pmaxT - failed to load package sprint")
        return()
    }

    if( !require("multtest", quietly=TRUE) ) {
        warning("Unit tests not run for pmaxT - failed to load package multtest")
        return()
    }

    # Find path for test scripts
    path <- system.file(package="sprint", "unitTests", "pmaxT")

    # Create test suite and execute it
    testSuite <- defineTestSuite(name="pmaxT", dirs=path)
    testData  <- runTestSuite(testSuite)
  
    # Print out the results
    printTextProtocol(testData, showDetails=TRUE)

#    tmp <- getErrors(testData)
#    if(tmp$nFail > 0 | tmp$nErr > 0) {
#        stop(paste("\n\npmaxT unit testing failed (#test failures: ", tmp$nFail,
#                   ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
#    }

}


sprint.test_pcor <- function() {

    if( !require("RUnit", quietly=TRUE) ) { 
        warning("Unit tests not run for pcor - failed to load package RUnit")
        return()
    }

    if( !require("sprint", quietly=TRUE) ) {
        warning("Uunit tests not run for pcor - failed to load package sprint")
        return()
    }

    # Find path for test scripts
    path <- system.file(package="sprint", "unitTests", "pcor")

    # Create test suite and execute it
    testSuite <- defineTestSuite(name="pcor", dirs=path)
    testData  <- runTestSuite(testSuite)
  
    # Print out the results
    printTextProtocol(testData, showDetails=TRUE)

}

