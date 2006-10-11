/******************************************************************************
Title: solar2mloci.cpp
Description: solar2mloci receives a string representing the name of a
             directory, a string representation of all the files in that
             directory, and the amount of those files as input.  createMloci
             converts all of the files in that directory to a sinlge mloci.out.
             This mloci.out file is a concatenation of the loci.out's created
             from the (m)ibds in the specified directory.  The individual
             loci.out's are separated in mloci.out by a line that contains a
             '#' in the first column and the name of the original (m)ibd file
             that the next information came from.  mloci.out is ouput in the
             (m)ibd directory specified in the first parameter.
Author: Eric Lunde, 10-23-03
Updates: Eric Lunde, 06-07-2005, I added a parameter, Sint *gzipWhenComplete.
         There is a parameter to delete.fixed.dir in solar2mloci.q.  If this is
         TRUE, there is no need to gzip the ibd files after the've been made
         into mloci.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cstring>
#include "verS.h"
using namespace std;

extern "C" {

void solar2mloci(char **directoryName, char **ibdNames, Sint *ibdNamesSize,
	       char **uniqueIds, Sint *gzipWhenComplete) {
  // Define the "gloabl" constants for the function
  const char *MLOCI_OUT = "mloci.out";
  char *MLOCI_OUT_GZ = (char *) malloc(sizeof(char) * strlen(MLOCI_OUT) + 4);
  strcpy(MLOCI_OUT_GZ, MLOCI_OUT);
  strcat(MLOCI_OUT_GZ, ".gz");
  //  const int FIELD_WIDTH = 15;
  const int END_OF_LINE_LENGTH = 8;
  const int COMMAND_LENGTH = 1024;
  const int PATH_LENGTH = 1024;
  
  // mloci_out is the handle to the output file MLOCI_OUT
  ofstream *mloci_out;

  // command is a char array to hold a system command
  char command[COMMAND_LENGTH];
  // ibdFile is the file we will be converting and appending to
  //   mloci.out.  Although the name of the variable is ibdFile, it
  //   could reference an mibd file as well
  ifstream *ibdFile;

  // The *Column variables will hold the information extracted from the
  // (m)ibd files
  int firstColumn, secondColumn;
  double thirdColumn, fourthColumn;

  // restOfLine will absorb the rest of the line after extracting
  //   the numerical values from the (m)ibd files
  char restOfLine[END_OF_LINE_LENGTH];

  // fullPathName will hold the directoryName concatenated with
  //   the file name.
  char fullPathName[PATH_LENGTH];

  // extension is a string of the last 3 characters of a file name to
  //   determine its extension type
  char *extension;

  // ibdFileCount is a bookkeeping variable to report to the user
  //   how many (m)ibd files created the MLOCI_OUT file
  int ibdFileCount = 0;

  // Construct the full path name for 'mloci.out'
  snprintf(fullPathName, PATH_LENGTH, "%s/%s", directoryName[0], MLOCI_OUT);
  mloci_out = new ofstream(fullPathName);

  // Test for successful opening
  if(mloci_out->fail()) {
    PROBLEM "The file %s could not be opened for writing.\nsolar2mloci.cpp key 105\n",
      fullPathName RECOVER(NULL_ENTRY);
  }

  cout << "Creating mloci.out:";

  // For each (non-directory) file in the directory, produce the loci.out
  // information and append it to mloci.out
  //  while( (dirEntry = readdir(directory)) != NULL ) {
  for(int i=0; i<*ibdNamesSize; i++) {
    if( (i % 35) == 0 ) {
      cout << endl;
    }
    cout << ". ";
    cout.flush();

    // Create the file's full (relative) path name to make referencing easier
    snprintf(fullPathName, PATH_LENGTH, "%s/%s", directoryName[0],
	     ibdNames[i]);

    snprintf(command, COMMAND_LENGTH, "test -d %s", fullPathName);

    // This if statement filters out directories
    if(system(command) == 0) {
      if(strcmp(ibdNames[i], ".") != 0
	 && strcmp(ibdNames[i], "..") != 0) {
	cout << "'" << ibdNames[i]  << "' is a directory -- ignored"
	     << endl;
      }
      continue;
    }

    // This if statement filters out mloci.out
    if(strcmp(ibdNames[i], MLOCI_OUT) == 0) {
      continue;
    }
    // This if statement filters out mloci.out.gz
    if(strcmp(ibdNames[i], MLOCI_OUT_GZ) == 0) {
      continue;
    }

    ibdFileCount++;

    // gunzip the (m)ibd file if necessary (this usually is necessary)
    extension = &fullPathName[ strlen(fullPathName)-3 ];
    if(strcmp(extension, ".gz") == 0) {
      snprintf(command, COMMAND_LENGTH, "gunzip %s", fullPathName);
      system(command);
      // Since the file is now gunzip'ed, we need to remove the '.gz' from
      // our file name
      fullPathName[ strlen(fullPathName)-3 ] = '\0';
      ibdNames[i][ strlen(ibdNames[i])-3 ] = '\0';
    }
    
    ibdFile = new ifstream(fullPathName);
    
    // Test for successful opening
    if(ibdFile->fail()) {
      PROBLEM "The file %s could not be opened for reading.\nsolar2mloci.cpp key 129\n",
	fullPathName RECOVER(NULL_ENTRY);
    }
    
    // Begin each new loci-appending with the '#' to indicate separation
    // and the name of the original file that the information came from
    *mloci_out << "# " << ibdNames[i] << endl;
    
    // Append to mloci.out the appropriate values.  If the first two columns
    // are NOT equal, the write the third column to mloci.out three times
    while(!ibdFile->eof()) {
      *ibdFile >> firstColumn >> secondColumn >> thirdColumn >> fourthColumn;
      if(firstColumn != secondColumn) {
	/*
	*mloci_out << setw(FIELD_WIDTH) << firstColumn
		   << setw(FIELD_WIDTH) << secondColumn
		   << setw(FIELD_WIDTH) << thirdColumn << endl;
		   */
	/*
	*mloci_out << firstColumn
		   << '\t' << secondColumn
		   << '\t' << thirdColumn << endl;
	*/
	*mloci_out << uniqueIds[firstColumn - 1]
		   << '\t' << uniqueIds[secondColumn - 1]
		   << '\t' << thirdColumn << endl;
	
      }
      ibdFile->getline(restOfLine, END_OF_LINE_LENGTH);
    }
    
    // Clean up the file references
    ibdFile->close();
    delete ibdFile;

    if(*gzipWhenComplete) {
      // Return the (m)ibd file to its gzip'ed form
      snprintf(command, COMMAND_LENGTH, "gzip -f %s", fullPathName);
      system(command);
    }
  }

  cout << "Done!" << endl;

  // Close and mloci.out
  delete mloci_out;

  snprintf(command, COMMAND_LENGTH, "mv -f %s/%s .", directoryName[0],
	   MLOCI_OUT);
  system(command);

  free(MLOCI_OUT_GZ);
}

/******************************************************************************
Title: alignMidbFiles
Description: alignMibdFiles is a preperation phase for createMloci.  In
             alignMibdFiles, we receive an array of file names that will be
             converted into mloci.out.  It extracts the last numerical value
             from mibd files and appends to the beginning of that numerical
             value "0"s to have all of the numerical values the same length
             while retaining their original magnitude.  This way createMloci
             will operate on them in alphabetcal order.
Input: char **directoryName - directoryName is the singular name of the
                              directory from which the (m)ibd files are
                              located.  Only the first index of this variable
                              is valid.
       char **ibdNames - ibdNames is the names of all the files currently in
                         the directory specified in directoryName.  (This may
                         not always be true, the user of this function could
                         limit the file names passed.  It really only needs to
                         hold one of the mibd file that has the most digits in
                         the last section of the file name and all of the other
                         file names with less digits in that same section.
       Sint *ibdNamesSize - ibdNamesSize is a scalar value representing the
                            amount of items in ibdNames.
Output: NONE
Side Effects: The file names specifed in ibdNames that have less then the most
              digits in the last section of the name will be altered on the
              disk.  The new file name (generated by this function) will be
              applied to those files.
Author: Eric Lunde, 10-24-03
******************************************************************************/
void alignMibdFiles(char **directoryName, char **ibdNames, Sint *ibdNamesSize)
{
  const int MAX_FILE_NAME_LENGTH = 128;
  const int MAX_DISTANCE_VALUE_LENGTH = 32;
  const int MAX_COMMAND_LENGTH = 256;

  int maxValueLength = -1;
  int amountOfZerosToAdd;
  char *prefix;
  char *chromNumber;
  char *distanceValue;
  char *restOfName;
  char tmpIbdName[MAX_FILE_NAME_LENGTH];
  char tmpDistanceValue[MAX_DISTANCE_VALUE_LENGTH];
  char command[MAX_COMMAND_LENGTH];

  // Get rid of potential trailing '/'
  if(directoryName[0][strlen(directoryName[0]) - 1] == '/') {
    directoryName[0][strlen(directoryName[0]) - 1] = '\0';
  }

  // Determine maximum length of the last section of the file name
  for(int i=0; i<*ibdNamesSize; i++) {
    // Copy the name of the file, becuase strtok will insert '\0's
    strncpy(tmpIbdName, ibdNames[i], MAX_FILE_NAME_LENGTH-1);
    snprintf(command, MAX_COMMAND_LENGTH, "test -d %s", tmpIbdName);

    // Extract the first token
    prefix = strtok(tmpIbdName, ".");

    // Test to see if the beginning of the file (pre ".") is "mibd" and is not
    // a directory
    if(system(command) != 0 && strcmp(prefix, "mibd") == 0) {
      // Extract the chromosome number and the distance value (the value that
      // will be used to determine if appending of "0"s will be necessary
      chromNumber = strtok(NULL, ".");
      distanceValue = strtok(NULL, ".");

      // Update the maximum amount of digits encountered if necessary
      if((int)strlen(distanceValue) > maxValueLength) {
	maxValueLength = strlen(distanceValue);
      }      
    }
  }

  // Update file names as necessary
  for(int i=0; i<*ibdNamesSize; i++) {
    // Copy the name of the file, becuase strtok will insert '\0's
    strncpy(tmpIbdName, ibdNames[i], MAX_FILE_NAME_LENGTH-1);
    snprintf(command, MAX_COMMAND_LENGTH, "test -d %s", tmpIbdName);

    // Extract the first token
    prefix = strtok(tmpIbdName, ".");

    // Test to see if the beginning of the file (pre ".") is "mibd" and is not
    // a directory
    if(system(command) != 0 && strcmp(prefix, "mibd") == 0) {
      // Extract the chromosome number, the distance value, and the rest of the
      // name
      chromNumber = strtok(NULL, ".");
      distanceValue = strtok(NULL, ".");
      restOfName = strtok(NULL, "\0");
      amountOfZerosToAdd = maxValueLength - (int)strlen(distanceValue);

      // If the file needs "0"s appended to the beginning, do so
      if(amountOfZerosToAdd > 0 ) {
	for(int j=0; j<amountOfZerosToAdd; j++) {
	  tmpDistanceValue[j] = '0';
	}
	tmpDistanceValue[ strlen(tmpDistanceValue) ] = '\0';
	strncat(tmpDistanceValue, distanceValue,
		MAX_DISTANCE_VALUE_LENGTH-1-strlen(distanceValue));

	// Generate and execute the system command to rename the file
	snprintf(command, MAX_COMMAND_LENGTH, "mv %s/%s %s/%s.%s.%s.%s",
		directoryName[0], ibdNames[i], directoryName[0], prefix,
		chromNumber, tmpDistanceValue, restOfName);
	system(command);
      }
    }  
  }
}

}
