/************************************************************************
 ************************************************************************
  fixibd.c
  Copyright (c) 2002 by MAYO CLINIC
  ALL RIGHTS RESERVED
************************************************************************
************************************************************************
 Date: February 18th, 2002                Coded by: Md. Asraful Isam
 Filename: fixibd.c
 Source file: fixibd.c
************************************************************************
 Input : input directory, output directory, name of phi file
 Output: corrected m/ibd files (that came from the input directory)
         in gzip format in the output directory

 Program description: 
 The program is to fix the ibd/mibd file with the help of phi2 file.
 Preprocess:
	- the fix all the files that is in m/ibd dir

 The steps to fix are: 
      1. Make the phi2 file ready (sort it with respect to the 2nd column) 
      2. Unzip all the files that is in the INPUT folder
      3. Read each file in the INPUT folder, 
		a. sort the file according to the column 2
		b. fix the file 
      4. Zip INPUT & OUTPUT dir 
	   

********************************************************************
 Development environment:
 Unix
********************************************************************
 Start date:February 18th, 2002 
 Update history:
           DATE                       MODIFICATIONS
         February 18, 2002            develop: get_phi_file_ready()
					develop: fix_mibd_file()
					develop: sort_mibd_file()
	 February 19, 2002		develop: main()
					develop: un/zip_dir()
6-17-03, Eric Lunde, I updated the error messages printed when files do not
         open.  I added line numbers to aid debugging and printed the messages
         to stderr instead of stdout
6-18-03, Eric Lunde, I added code in the fix_ibd_file function to test the
         existance of the output directory before files were to be placed
         there.  If the directory did not exist, it now creates that directory.
8-20-03, Eric Lunde, I added code in fixibd to disregard (do not fix)
         subdirectories.  The purpose was to test fixibd on only one file, with
         all the others in a subdirectory without the program crashing trying
         to access the subdirectory.
9-16-03, Eric Lunde, I got rid of the global variables.  Now everything is
         local.  This in turn required redefining the parameter lists for the
         utility functions.  I added code to get rid of an extension of '.gz'
         for the phi file if the '.gz' didn't exist.  I also updated the
         function documentation.

 Test protocol/history:
 Manually tested & run successfully
********************************************************************
 Note:
 - The current version of Solar creates IBD matrices without zeros 
   for missing data.  This program reads in the MIBD/IBD files 
   and adds the zeros back in.
 - make sure the out out dir (e.g. mibd_new) is created
          
********************************************************************
********************************************************************/
/* Library */
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <iostream>
#include "verS.h"
using namespace std;
 
/* Define */
#define DIRECTORY_NAME_LENGTH 128
#define FILE_NAME_LENGTH 128
#define COMMAND_LENGTH 512
#define LINE_LENGTH 1024

/* Forward declarations for functions */
void create_sorted_phi2_file(char *phiFileName, char *sortedPhiFileName);
void unzip_dir(char *directoryName);
void zip_dir(char *directoryName);
void sort_ibd_file(char *ibdDirectoryName, char *file_name,
		   char *temporaryIbdName);
void fix_ibd_file(char *ibdDirectoryName,  char *ibd_file,
		  char *fixedIbdDirectoryName, char *sortedPhiFileName,
		  char *temp_ibd);

extern "C" {

/*****************************
File name: fixibd.c
Method name: fixibd
Description: (see the desctiption at the top of this file)
Input: char **ibd_directory_name - This 2 dimensional array of chars is only
                                   used for the first array's value.  This
                                   string specifies the name of the folder to
                                   be fix'ed.  This value can be relative to
                                   the executing file or absolute from the
                                   root.
       char **fixed_ibd_directory_name - This parameter acts similarly to
                                         ibd_directory_name except that this
                                         string value described the name of the
                                         output directory where the fix'ed
                                         ibd's will be placed.
       char **phi_file_name - Similar in format to the other two variables,
                              phi_file_name specifies the name (and location)
                              of the phi file used to fix the ibds.
Output: NONE
Side Effects: If the output directory does not exist, this function creates it.
Author: Eric Lunde, 7-08-03
******************************/
void fixibd(char **ibd_directory_name, char **fixed_ibd_directory_name,
	    char **phi_file_name, char **ibd_directory_name_to_print)
{
  char ibdDirectoryName[DIRECTORY_NAME_LENGTH];
  char fixedIbdDirectoryName[DIRECTORY_NAME_LENGTH];
  char phiFileName[FILE_NAME_LENGTH];
  char sortedPhiFileName[FILE_NAME_LENGTH];
  char temporaryIbdName[FILE_NAME_LENGTH];
  char *ibdDirectoryNameToPrint
    = (char *) malloc(strlen(*ibd_directory_name_to_print) + 1);
  strcpy(ibdDirectoryNameToPrint, *ibd_directory_name_to_print);

  char command[COMMAND_LENGTH];
  char *extension;
  int ibdFileCounter = 0;

  /* Directory interface variables */
  struct  dirent  *dptr;
  DIR   *dirp;

  /* Copy file and directory names to local memory */
  strncpy(ibdDirectoryName, ibd_directory_name[0], DIRECTORY_NAME_LENGTH-1);
  strncpy(fixedIbdDirectoryName, fixed_ibd_directory_name[0],
	  DIRECTORY_NAME_LENGTH-1);
  strncpy(phiFileName, phi_file_name[0], FILE_NAME_LENGTH-1);

  /* Remove the trailing '/' (if it exists) from the directory names for
     consistancy */
  if(ibdDirectoryName[strlen(ibdDirectoryName)-1] == '/') {
    ibdDirectoryName[strlen(ibdDirectoryName)-1] = '\0';
  }
  if(fixedIbdDirectoryName[strlen(fixedIbdDirectoryName)-1] == '/') {
    fixedIbdDirectoryName[strlen(fixedIbdDirectoryName)-1] = '\0';
  }

  /* Test if the phi file is gzip'ed */
  extension = &phiFileName[strlen(phiFileName)-3];  
  if(strcmp(extension, ".gz") == 0) {
    snprintf(command, COMMAND_LENGTH-1, "test -r %s", phiFileName);
    /* If the file exists and is readable */
    if(system(command) == 0) {
      snprintf(command, COMMAND_LENGTH-1, "gunzip -f %s", phiFileName);
      system(command);
    }
    *extension = '\0';
  }

  /* Generate sorted phi file name and create it*/
  strncpy(sortedPhiFileName, phiFileName, FILE_NAME_LENGTH-1);
  strncat(sortedPhiFileName, ".sorted",
	  FILE_NAME_LENGTH-strlen(sortedPhiFileName)-1);
  create_sorted_phi2_file(phiFileName, sortedPhiFileName);

  /* Generate the temporary ibd file name */
  strncpy(temporaryIbdName, ibdDirectoryName, FILE_NAME_LENGTH-1);
  /* 10 is the length of the appended string */
  strncat(temporaryIbdName, "/temp_ibd",
	  FILE_NAME_LENGTH-10-1);

  /* unzip all files of the input dir */
  printf("Unzipping files in the input dir: %s\n", ibdDirectoryNameToPrint);
  free(ibdDirectoryNameToPrint);
  unzip_dir(ibdDirectoryName);
  
  /* get access to directory structure */
  if((dirp=opendir(ibdDirectoryName))==NULL) { 
    PROBLEM "Error opening %s\nfixibd.c key 116\n",
      ibdDirectoryName RECOVER(NULL_ENTRY);
  }
  
  /* read each file in the directory */
  printf("Fixing files:");
  while( (dptr=readdir(dirp)) ) {
    if( strcmp(dptr->d_name, "temp_ibd") == 0) {
      continue;
    }

    sprintf(command, "test -d %s/%s", ibdDirectoryName, dptr->d_name);
    /* If the file is not a directory */
    if(system(command) == 0) {
      continue;
    }

    if(strncmp(dptr->d_name, "mibd.", strlen("mibd.")) != 0
       && strncmp(dptr->d_name, "ibd.", strlen("ibd.")) != 0) {
      continue;
    }

    if( (ibdFileCounter % 35) == 0 ) {
      printf("\n");
    }
    printf(". ");
    fflush(stdout);
    ibdFileCounter++;
    
    sort_ibd_file(ibdDirectoryName, dptr->d_name, temporaryIbdName);
    fix_ibd_file(ibdDirectoryName, dptr->d_name, fixedIbdDirectoryName,
		 sortedPhiFileName, temporaryIbdName);
  }
  closedir(dirp);
  printf("Done!\n");
  

  sprintf(command, "\\rm -f %s %s", temporaryIbdName, sortedPhiFileName);
  system(command);

  /* zip all files of the given input directory */
  /* This no longer has to be done because we are copying all (m)ibd files
     to a local copy to operate.  They will be deleted at the end, so there
     is no need to waste time to gzip them when we are going to delete them
     very very soon.  Eric Lunde 2005/07/27
     printf("Zipping files in the input dir: %s\n", ibdDirectoryName);
     zip_dir(ibdDirectoryName); */
    

  /* zip all files of the given output directory */
  /*
    printf("Zipping files in the output dir: %s", fixedIbdDirectoryName);
    zip_dir(fixedIbdDirectoryName);
    printf("...Done!\n");
    */
}

}

/************************************************************************
 Function: create_sorted_phi2_file
 ************************************************************************
 Input : The phi file (which is made up of the path) and the name for the
         sorted phi file.
 Output: NONE
 Side Effects: 'phiFileName'.sorted is created and is a sorted version of
               'phiFileName' based on the second column of the file

 Function description:
 This function is to sort the phi file with respect to the column #2
 the sorted phi file is going to be used to fix the m/ibd files 
************************************************************************/
void create_sorted_phi2_file(char *phiFileName, char *sortedPhiFileName)
{
  char sort_command[COMMAND_LENGTH];

  // Adding the -n switch should fix a bug on 64-bit systems.
  // Also now using -k flag instead of old way of specifying keys.
  sprintf(sort_command, "sort -k 2,2n -k 1,1n %s -o %s", phiFileName,
	  sortedPhiFileName);
  system(sort_command);
}

/************************************************************************
 Function: unzip_dir
************************************************************************
 Input : dir to unzip all of the files in it 
 Output: unzip all the files in a given dir
 Return: none

 Function description:
 to unzip all of the files in a given dir, using shell command gunzip 
************************************************************************/
void unzip_dir(char *dir_name)
{
  char command[COMMAND_LENGTH];
  
  sprintf(command, "test -d %s", dir_name);
  // If the system command exits with failure (non zero, or true in C++),
  // issue an error message and quit.
  if( system(command) ) {
    PROBLEM "The directory '%s' is not a directory.\nfixibd.cpp key 262\n",
      dir_name RECOVER(NULL_ENTRY);
  }

  /* execute shell command */
  sprintf(command,"gunzip -qf %s/*", dir_name);
  int systemResult = system(command);
  if( systemResult ) {
    cout << "No files were unzipped in directory: " << dir_name << endl;
  }
}

/***********************************************************************
 Function: zip_dir
************************************************************************
 Input : dir to zip all of the files in it 
 Output: zip all the files in a given dir
 Return: none

 Function description:
 to zip all of the files in a given dir, using shell command gzip 
************************************************************************/
void zip_dir(char *dir_name)
{
  char zip_command[COMMAND_LENGTH];
  
  sprintf(zip_command,"gzip -f %s/*.*", dir_name);

  /* executing shell command */
  system(zip_command);
}

/************************************************************************
 Function: sort_mibd_file
************************************************************************
 Input : ibd directory, file name, and temporary file name
 Output: sorted m/ibd file with the name of temp_ibd
 Return: None

 Function description: 
 This function is to sort the m/imd file with respect to the column #2
************************************************************************/
void sort_ibd_file(char *ibdDirectoryName, char *file_name,
		   char *temporaryIbdName)
{
  char sort_command[COMMAND_LENGTH];
  char file_path[FILE_NAME_LENGTH];
  
  sprintf(file_path, "%s/%s", ibdDirectoryName, file_name);

  // Adding the -n switch should fix a bug on 64-bit systems
  // Also now using -k flag instead of old way of specifying keys.
  sprintf(sort_command,"sort -k 2,2n -k 1,1n %s -o %s", file_path, temporaryIbdName);
  system(sort_command);
}

/************************************************************************
 Function: fix_ibd_file
************************************************************************
 Input : m/ibd file (comprised of the ibd directory and the specific file name)
         to fix, the destination (fixed) directory, sorted phi file, and the
         temporary file name
 Output: fixed m/ibd file in the given output dir 
 Return: None

 Function description: 
 Function to fix the m/ibd files using phi file
 NOTE:
 - m/ibd files are incomplete with some missing rows
 - missing rows are present in phi file
 - keep exixting rows in the m/ibd files
 - find rows that are missing in m/ibd
 - merge then together   
************************************************************************/
void fix_ibd_file(char *ibdDirectoryName,  char *ibd_file,
		  char *fixedIbdDirectoryName, char *sortedPhiFileName,
		  char *temp_ibd)
{
  FILE* fp_phi;
  FILE* fp_mibd;
  FILE* fp_new_mibd;
  
  int ph1,ph2,mi1,mi2, previousMi1, previousMi2;
  double ph3,ph4,mi3,mi4;
  int flag=0;
  char output_file[FILE_NAME_LENGTH];
  char command[COMMAND_LENGTH];
  char line[LINE_LENGTH];
  char line2[LINE_LENGTH];
  
  fp_mibd = fopen(temp_ibd, "r");
  if(!fp_mibd) {
    PROBLEM "Error opening file: %s.\nfixibd.c key 260\n",
      temp_ibd RECOVER(NULL_ENTRY);
  }
  
  fp_phi = fopen(sortedPhiFileName, "r");
  if(!fp_phi) {
    PROBLEM "Error opening file: %s.\nfixibd.c key 267\n",
      sortedPhiFileName RECOVER(NULL_ENTRY);
  }

  /**********************
  This next sprintf and if branch were added to create the output directory
  (if it didn't already exist) before the program writes to its contents.
  The system command inside of the if conditional returns 0 if the file exists
  and is a directory, a non-zero answer otherwise.
  Eric Lunde, 6-18-03
  **********************/
  sprintf(command, "test -d %s", fixedIbdDirectoryName);  
  if(system(command) != 0 ) {
    sprintf(command, "mkdir %s", fixedIbdDirectoryName);
    system(command);
  }
  
  sprintf(output_file,"%s/%s",fixedIbdDirectoryName, ibd_file);
  fp_new_mibd = fopen(output_file,"w");
  if(!fp_new_mibd) {
    PROBLEM "Error opening file: %s.\nMake an output dir: %s.\nfixibd.c key 275\n",
      output_file, fixedIbdDirectoryName RECOVER(NULL_ENTRY);
  }

  previousMi1 = previousMi2 = -1;
  fgets(line, LINE_LENGTH-1, fp_mibd);
  while(sscanf(line,"%d %d %lf %lf",&mi1,&mi2,&mi3,&mi4) == 4
	&& !feof(fp_mibd)) {
    flag=0;
    if(previousMi1 == mi1 && previousMi2 == mi2) {
      cerr << endl << "Warning!  The ids: " << mi1 << " and " << mi2
	   << " were found in consecutive rows." << endl
	   << "This is not a valid (m)ibd file.  Processing will continue."
	   << endl << endl;
    }

    fgets(line2, LINE_LENGTH-1, fp_phi);
    int scanResult = sscanf(line2,"%d %d %lf %lf",&ph1,&ph2,&ph3,&ph4);
    //    while(fscanf(fp_phi,"%d %d %lf %lf",&ph1,&ph2,&ph3,&ph4) == 4) {
    while(scanResult == 4 && !feof(fp_phi)) {
      if((mi2==ph2) && (mi1==ph1)) {
	fprintf(fp_new_mibd,"%4d %4d %f %f\n", mi2, mi1, mi3, mi4);
	flag=1;
      }else {
	fprintf(fp_new_mibd,"%4d %4d %f %f\n", ph2, ph1, ph3, ph4);
      }
      if(flag==1) {
	break;
      }
      fgets(line2, LINE_LENGTH-1, fp_phi);
      scanResult = sscanf(line2,"%d %d %lf %lf",&ph1,&ph2,&ph3,&ph4);
    }
    if(flag==0) {
      break;
    }

    previousMi1 = mi1;
    previousMi2 = mi2;
    fgets(line, LINE_LENGTH-1, fp_mibd);
  }
  /*Closing all the opened files*/
  fclose(fp_phi);
  fclose(fp_mibd);
  fclose(fp_new_mibd);
  return;
}
