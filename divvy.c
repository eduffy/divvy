
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <pcre.h>

#define MAX_PATH_LENGTH  (4096)
#define MAX_PARTIAL_SIZE (1024)
#define TAG_PARTIAL      (0)

#define streq(X,Y)      (strcmp((X),(Y))==0)
#define startswith(X,P) (strncmp((X), (P), strlen(P)) == 0)


typedef struct {
  char *name;
  char *pattern;
} predefine_t;

predefine_t PREDEFINE_LIST[] =
{
   "--fastq", "^@.*\\n.*\\n\\+",
   NULL,      NULL
};

struct buffer {
  char *data, *start, *end;
};

int file_exists(const char *path)
{
   return access(path, R_OK) == 0;
}

long getfilesize(const char *path)
{
  struct stat buf;
  if(stat(path, &buf) == -1)
    return -1;
  return buf.st_size;
}

int parse_commandline(int argc, char **argv,
                      char **pattern,
                      char **infile)
{
   while(--argc) {
      ++argv;
      if(streq(argv[0], "-r")) {
         *pattern = strdup(argv[1]);
         --argc;
      }
      else if(startswith(argv[0], "-r")) {
         *pattern = strdup(argv[0] + 2);
      }
      else if(startswith(argv[0], "--regex=")) {
         *pattern = strdup(argv[0] + 8);
      }
      else if(startswith(argv[0], "--")) {
         predefine_t *p = &PREDEFINE_LIST[0];
         for(; p->name != NULL; ++p) {
            if(streq(p->name, argv[0])) {
               *pattern = strdup(p->pattern);
            }
         }
      }
      else {
         if(file_exists(argv[0])) {
            *infile = strdup(argv[0]);
         }
      }
   }
}
void advance_record(const char *pattern, struct buffer *buf)
{
   pcre       *regex = NULL;
   const char *err;
   int         erroffset;
   int         rank;
   int         retcode;
   int         matches[3] = { 0 };

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   regex = pcre_compile(pattern, PCRE_MULTILINE, &err, &erroffset, NULL);
   if(regex == NULL) {
      fprintf(stderr, "An error occured while compiling regular expression, "
                      "%s at offset %d.\n", err, erroffset);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   
   retcode = pcre_exec(regex,
                       NULL,
                       buf->start,
                       buf->end - buf->start,
                       0, 0, matches, 3);
   if(retcode < 0) {
      fprintf(stderr, "No record header found in chunk %d "
                      "matching regular expression %s.\n", 
                      rank, pattern);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   buf->start += matches[0];
}

void read_fastq(char *filename, struct buffer *buf)
{
   MPI_File fh;
   MPI_Status status;
   int    size, rank;
   long   filesize, chunksize;
   int    count;
   int    rescode;

   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   filesize  = getfilesize(filename);
   chunksize = (long) ceil((double)filesize / size);

   if(rank == 0) {
      printf("File size:  %12ld bytes\n", filesize);
      printf("Chunk size: %12ld bytes\n", chunksize);
      printf("Num chunks: %12d\n", size);
   }

   rescode = MPI_File_open(MPI_COMM_WORLD,
                           filename,
                           MPI_MODE_RDONLY,
                           MPI_INFO_NULL,
                           &fh);
   if(rescode == MPI_ERR_IO) {
      fprintf(stderr,
              "Process #%d cannot open %s for reading.\n",
              rank, filename);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   buf->data = malloc(chunksize + MAX_PARTIAL_SIZE);
   MPI_File_set_view(fh,
                     rank * chunksize,
                     MPI_CHAR, MPI_CHAR,
                     "native", MPI_INFO_NULL);
   MPI_File_read(fh, buf->data, chunksize, MPI_CHAR, &status);
   MPI_Get_count(&status, MPI_CHAR, &count);
   MPI_File_close(&fh);

   buf->start = buf->data;
   buf->end   = buf->data + count;
}

void transfer_partials(const char *pattern, struct buffer *buf)
{
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if(rank > 0) {
      if(pattern != NULL)
         advance_record(pattern, buf);
      MPI_Send(buf->data,
               buf->start - buf->data,
               MPI_CHAR,
               rank - 1, TAG_PARTIAL, MPI_COMM_WORLD);
   }
   if(rank < size - 1) {
      MPI_Status stat;
      int count;

      MPI_Recv(buf->end, MAX_PARTIAL_SIZE, MPI_CHAR,
               rank + 1, TAG_PARTIAL, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_CHAR, &count);
      buf->end += count;
   }
}

void write_chunks(char *filename, struct buffer *buf)
{
   int rank;
   char outfilename[MAX_PATH_LENGTH];
   FILE *fh;

   /* FIXME: use MPI_File */

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   sprintf(outfilename, "%s.%03d", filename, rank);
   fh = fopen(outfilename, "w");
   fwrite(buf->start, 1, buf->end - buf->start, fh);
   fclose(fh);
}

int main(int argc, char *argv[])
{
   struct buffer  buf;
   char          *pattern = NULL;
   char          *infile  = NULL;

   parse_commandline(argc, argv, &pattern, &infile);
   if(infile == NULL)
      return EXIT_FAILURE;

   MPI_Init(&argc, &argv);

   read_fastq(infile, &buf);
   transfer_partials(pattern, &buf);
   write_chunks(infile, &buf);

   MPI_Finalize();
   free(buf.data);
   free(pattern);
   free(infile);
   return EXIT_SUCCESS;
}
