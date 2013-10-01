/*
 * Copyright (C) 2013 Edward B. Duffy
 * License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
 */

#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <glob.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

#define TAG_CHUNK       (0)

#define streq(X,Y)      (strcmp((X),(Y))==0)
#define startswith(X,P) (strncmp((X), (P), strlen(P)) == 0)

#define file_exists(P)  (access((P), R_OK) == 0)

long getfilesize(const char *path)
{
   struct stat buf;
   if(stat(path, &buf) == -1)
      return -1;
   return buf.st_size;
}

int parse_commandline(int argc, char **argv,
                      char **pattern,
                      char **outfile)
{
   while(--argc) {
      ++argv;
      if(streq(argv[0], "-p")) {
         *pattern = strdup(argv[1]);
         --argc;
         ++argv;
      }
      else if(streq(argv[0], "-o")) {
         *outfile = strdup(argv[1]);
         --argc;
         ++argv;
      }
/*
      else {
         if(file_exists(argv[0])) {
            *infile = strdup(argv[0]);
         }
      }
*/
   }
}

int main(int argc, char *argv[])
{
   MPI_File  fh;
   MPI_Status status;
   FILE     *fp;
   char     *buffer;
   char     *pattern;
   char     *outfile;
   glob_t    pglob;
   int       globrv;
   int       rank, size;
   long      filesize;
   long      offset;
   int       rescode;

   MPI_Init(&argc, &argv);
   parse_commandline(argc, argv, &pattern, &outfile);

   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   globrv = glob(pattern, 0, NULL, &pglob);
   if(rank == 0 && pglob.gl_pathc < size) {
     fprintf(stderr, "Too few files given; %d < %d\n", pglob.gl_pathc, size);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   else if(rank == 0 && pglob.gl_pathc > size) {
     fprintf(stderr, "Too many files given; %d > %d\n", pglob.gl_pathc, size);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }

   filesize = getfilesize(pglob.gl_pathv[rank]);
   buffer = malloc(filesize);
   assert(buffer != NULL);
   fp = fopen(pglob.gl_pathv[rank], "r");
   assert(fp != NULL);
   fread(buffer, filesize, 1, fp);
   fclose(fp);

   MPI_Scan(&filesize, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
   offset -= filesize;
   printf("rank = %03d; offset = %12ld\n", rank, offset);

   rescode = MPI_File_open(MPI_COMM_WORLD,
                           outfile,
                           MPI_MODE_CREATE | MPI_MODE_WRONLY,
                           MPI_INFO_NULL,
                           &fh);
   if(rescode == MPI_ERR_IO) {
      fprintf(stderr,
              "Process #%d cannot open %s for reading.\n",
              rank, outfile);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   MPI_File_set_view(fh,
                     offset,
                     MPI_CHAR, MPI_CHAR,
                     "native", MPI_INFO_NULL);
   MPI_File_write(fh, buffer, filesize, MPI_CHAR, &status);
   MPI_File_close(&fh);

   MPI_Finalize();
   return EXIT_SUCCESS;
}
