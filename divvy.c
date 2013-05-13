
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

#define MAX_PATH_LENGTH  (4096)
#define MAX_PARTIAL_SIZE (1024)
#define TAG_PARTIAL      (0)

struct buffer {
  char *data, *start, *end;
};

long getfilesize(const char *fn)
{
  struct stat buf;
  if(stat(fn, &buf) == -1)
    return -1;
  return buf.st_size;
}

void locate_next_sequence(struct buffer *buf)
{
   char *p;
   while(buf->start[0] != '@') ++(buf->start);
   p = buf->start + 1;
   while(p[0] != '\n' || p[1] != '+') {
      if(p[0] == '@')
         buf->start = p;
      ++p;
   }
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

void transfer_partials(struct buffer *buf)
{
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if(rank > 0) {
      locate_next_sequence(buf);
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
   struct buffer buf;

   MPI_Init(&argc, &argv);

   read_fastq(argv[1], &buf);
   transfer_partials(&buf);
   write_chunks(argv[1], &buf);

   MPI_Finalize();
   free(buf.data);
   return EXIT_SUCCESS;
}
