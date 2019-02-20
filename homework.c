#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
typedef unsigned char uchar;

float smooth[3][3] = {
  {1 / 9.0f, 1 / 9.0f, 1 / 9.0f},
  {1 / 9.0f, 1 / 9.0f, 1 / 9.0f},
  {1 / 9.0f, 1 / 9.0f, 1 / 9.0f}
};

float blur[3][3] = {
  {1 / 16.0f, 2 / 16.0f, 1 / 16.0f},
  {2 / 16.0f, 4 / 16.0f, 2 / 16.0f},
  {1 / 16.0f, 2 / 16.0f, 1 / 16.0f}
};

float sharpen[3][3] = {
  {0, -2 / 3.0f, 0},
  {-2 / 3.0f, 11 / 3.0f, -2 / 3.0f},
  {0, -2 / 3.0f, 0}
};

float mean[3][3] = {
  {-1.0f, -1.0f, -1.0f},
  {-1.0f, 9.0f, -1.0f},
  {-1.0f, -1.0f, -1.0f}
};

float emboss[3][3] = {
  {0.0f, 1.0f, 0.0f},
  {0.0f, 0.0f, 0.0f},
  {0.0f, -1.0f, 0.0f}
};

char *bssembssem[] = {"blur", "smooth", "sharpen",
                     "emboss", "mean", "blur", "smooth",
                     "sharpen", "emboss", "mean"};

typedef struct {
  uchar red;
  uchar green;
  uchar blue;
} pixel;

typedef struct {
  float red;
  float green;
  float blue;
} threeFloats;

typedef struct {
  char fileType[3];
  char width[6], height[6];
  char maxVal[4];
  uchar **grayScale;
  pixel **color;
} image;

void readInput(const char * fileName, image *img) {
  FILE *file = fopen(fileName, "rb");
  char read;
  int idx = 0, boolean = 0;

  // Set all to 0
  memset(img->fileType, 0, 3);
  memset(img->width, 0, 5);
  memset(img->height, 0, 5);
  memset(img->maxVal, 0, 4);

  // Read the "magic number"(P5 or P6)
  fread(img->fileType, 2, sizeof(char), file);
  fread(&read, sizeof(char), 1, file);
  fread(&read, sizeof(char), 1, file);

  // Read the width and height of the image
  while (read != '\n') {
    if (read == ' ') {
      boolean = 1;
      idx = 0;
      fread(&read, sizeof(char), 1, file);
      continue;
    }
    if (!boolean) {
      img->width[idx] = read;
    } else {
      img->height[idx] = read;
    }
    fread(&read, sizeof(char), 1, file);
    idx++;
  }

  int height = atoi(img->height);
  int width = atoi(img->width);
  idx = 0;
  fread(&read, sizeof(char), 1, file);

  // Read the maximum color value(img->maxVal)
  while (read != '\n') {
    img->maxVal[idx] = read;
    fread(&read, sizeof(char), 1, file);
    idx++;
  }

  // Found a .pgm image. Save the pixels in a matrix of char(img->grayScale)
  if (strstr(img->fileType, "P5")) {
    img->grayScale = (uchar **) malloc(height * sizeof(uchar *));
    img->color = NULL;
    for (int i = 0; i < height; i++) {
      img->grayScale[i] = (uchar *) malloc(width * sizeof(uchar));
      fread(img->grayScale[i], sizeof(uchar), width, file);
    }
  }

  // Found a .pnm image. Save the pixels in a matrix(img->color)
  else {
    img->grayScale = NULL;
    img->color = (pixel **) malloc(height * sizeof(pixel *));
    for (int i = 0; i < height; i++) {
      img->color[i] = (pixel *) malloc(width * sizeof(pixel));
      fread(img->color[i], sizeof(pixel), width, file);
    }
  }
  fclose(file);
}

void writeData(const char * fileName, image *img) {
  FILE *file = fopen(fileName, "wb");
  char newLine = '\n', space = ' ';

  // Write the header of the image
  fwrite(img->fileType, 1, strlen(img->fileType), file);
  fwrite(&newLine, 1, 1, file);
  fwrite(img->width, 1, strlen(img->width), file);
  fwrite(&space, 1, 1, file);
  fwrite(img->height, 1, strlen(img->height), file);
  fwrite(&newLine, 1, 1, file);
  fwrite(img->maxVal, 1, strlen(img->maxVal), file);
  fwrite(&newLine, 1, 1, file);

  int height = atoi(img->height);
  int width = atoi(img->width);

  // Write the pixels of the image based on its type (P5 or P6)
  if (strstr(img->fileType, "P5")) {
    for (int i = 0; i < height; i++) {
      fwrite(img->grayScale[i], 1, width, file);
      free(img->grayScale[i]);
    }
    free(img->grayScale);
  }
  else {
    for (int i = 0; i < height; i++) {
      fwrite(img->color[i], 3, width, file);
      free(img->color[i]);
    }
    free(img->color);
  }
  fclose(file);
}

void applyFilter(void ***matrix, const char *fileType,
  const char *filterName, int height, int width) {

  float K[3][3];
  if (filterName == NULL)
    return;

  // Assign K to the corresponding filter
  if (!strcmp(filterName, "smooth")) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        K[i][j] = smooth[i][j];
  } else if (!strcmp(filterName, "blur")) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        K[i][j] = blur[i][j];
  } else if (!strcmp(filterName, "sharpen")) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        K[i][j] = sharpen[i][j];
  } else if (!strcmp(filterName, "mean")) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        K[i][j] = mean[i][j];
  } else if (!strcmp(filterName, "emboss")) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        K[i][j] = emboss[i][j];
  } else {
    return;
  }


  if (!strcmp(fileType, "P5")) {
    uchar **grayScale = (uchar **) *matrix;
    uchar **aux = (uchar **) malloc(sizeof(uchar *) * height);
    for (int i = 0; i < height; ++i)
      aux[i] = (uchar *) malloc (sizeof(uchar) * width);

    // Make a copy of the original matrix
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        aux[i][j] = grayScale[i][j];
      }
    }

    for (int i = 1; i < height - 1; ++i) {
      for (int j = 1; j < width - 1; ++j) {
        float result[3][3];

        // Multiply each element of the 2 matrixes
        for (int k = i - 1; k < i + 2; ++k) {
          for (int l = j - 1; l < j + 2; ++l) {
            result[k - (i - 1)][l - (j - 1)] = aux[k][l]
              * K[k - (i - 1)][l - (j - 1)];
          }
        }

        // Add all elements from the result matrix
        float filter = 0;
        for (int k = 0; k < 3; ++k) {
          for (int l = 0; l < 3; ++l) {
            filter += result[k][l];
          }
        }
        grayScale[i][j] = filter;
      }
    }

    // Free the copy
    for (int i = 0; i < height; ++i)
      free(aux[i]);
    free(aux);

  } else {
    pixel **color = (pixel **) *matrix;
    pixel **aux = (pixel **) malloc(sizeof(pixel *) * height);
    for (int i = 0; i < height; ++i)
      aux[i] = (pixel *) malloc(sizeof(pixel) * width);

    // Make a copy of the original matrix
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        aux[i][j].red = color[i][j].red;
        aux[i][j].green = color[i][j].green;
        aux[i][j].blue = color[i][j].blue;
      }
    }

    for (int i = 1; i < height - 1; ++i) {
      for (int j = 1; j < width - 1; ++j) {
        threeFloats result[3][3];

        // Multiply each element of the 2 matrixes
        for (int k = i - 1; k < i + 2; ++k) {
          for (int l = j - 1; l < j + 2; ++l) {
            result[k - (i - 1)][l - (j - 1)].red = aux[k][l].red
              * K[k - (i - 1)][l - (j - 1)];
            result[k - (i - 1)][l - (j - 1)].green = aux[k][l].green
              * K[k - (i - 1)][l - (j - 1)];
            result[k - (i - 1)][l - (j - 1)].blue = aux[k][l].blue
              * K[k - (i - 1)][l - (j - 1)];
          }
        }

        // Add all elements from the result matrix
        threeFloats filter;
        filter.red = 0;
        filter.green = 0;
        filter.blue = 0;
        for (int k = 0; k < 3; ++k) {
          for (int l = 0; l < 3; ++l) {
            filter.red += result[k][l].red;
            filter.green += result[k][l].green;
            filter.blue += result[k][l].blue;
          }
        }
        color[i][j].red = filter.red;
        color[i][j].green = filter.green;
        color[i][j].blue = filter.blue;
      }
    }

    // Free the copy
    for (int i = 0; i < height; ++i)
      free(aux[i]);
    free(aux);
  }
}

void operationProc0(image *input, int nProcesses,
  char *filter, MPI_Datatype MPI_PIXEL) {

  int heightImage = atoi(input->height);
  int height = heightImage / 3;
  height *= 3;
  height = height / nProcesses;

  int heightOneProc = height;
  if (heightOneProc + 2 >= heightImage)
    heightOneProc = heightImage - 2;

  height = atoi(input->height);
  int width = atoi(input->width);

  if (!strcmp(input->fileType, "P6")) {
    // Send the submatrixes to the slave proccesses
    for (int i = 1; i < nProcesses; ++i) {
      if (i != nProcesses - 1)
        for (int j = i * heightOneProc; j < (i + 1) * heightOneProc + 2; ++j) {
          MPI_Send(input->color[j], width, MPI_PIXEL, i, 0, MPI_COMM_WORLD);
        }
      else {
        for (int j = i * heightOneProc; j < height; ++j) {
          MPI_Send(input->color[j], width, MPI_PIXEL, i, 0, MPI_COMM_WORLD);
        }
      }
    }

  } else if (!strcmp(input->fileType, "P5")) {
    // Send the submatrixes to the slave proccesses
    for (int i = 1; i < nProcesses; ++i) {
      if (i != nProcesses - 1)
        for (int j = i * heightOneProc; j < (i + 1) * heightOneProc + 2; ++j) {
          MPI_Send(input->grayScale[j], width, MPI_UNSIGNED_CHAR,
            i, 0, MPI_COMM_WORLD);
        }
      else {
        for (int j = i * heightOneProc; j < height; ++j) {
          MPI_Send(input->grayScale[j], width, MPI_UNSIGNED_CHAR,
            i, 0, MPI_COMM_WORLD);
        }
      }
    }
  }

  // Apply the filter on the proccess submatrix
  if (!strcmp(input->fileType, "P5"))
    applyFilter((void ***) &input->grayScale, input->fileType,
      filter, heightOneProc + 2, width);
  else if (!strcmp(input->fileType, "P6"))
    applyFilter((void ***) &input->color, input->fileType,
      filter, heightOneProc + 2, width);

  // Receive the results from the rest of the proccesses
  for (int i = 1; i < nProcesses; ++i) {
    if (i != nProcesses - 1) {
      if (!strcmp(input->fileType, "P5")) {
        int start = i * heightOneProc + 1;
        int end = (i + 1) * heightOneProc + 1;
        for (int j = start; j < end; ++j) {
          MPI_Recv(input->grayScale[j], width, MPI_UNSIGNED_CHAR,
            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      } else if (!strcmp(input->fileType, "P6")) {
        int start = i * heightOneProc + 1;
        int end = (i + 1) * heightOneProc + 1;
        for (int j = start; j < end; ++j) {
          MPI_Recv(input->color[j], width, MPI_PIXEL,
            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }

    } else {
      if (!strcmp(input->fileType, "P5")) {
        for (int j = i * heightOneProc + 1; j < height - 1; ++j) {
          MPI_Recv(input->grayScale[j], width, MPI_UNSIGNED_CHAR,
            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      } else if (!strcmp(input->fileType, "P6")) {
        for (int j = i * heightOneProc + 1; j < height - 1; ++j) {
          MPI_Recv(input->color[j], width, MPI_PIXEL,
            i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
    }
  }
}

void operationSlaveProc(int height, int width,
  char *fileType, char *filter, MPI_Datatype MPI_PIXEL) {

  uchar **grayScale;
  pixel **color;
  if (!strcmp(fileType, "P5")) {
    color = NULL;
    grayScale = (uchar **) malloc(sizeof(uchar *) * height);
    for (int i = 0; i < height; ++i)
      grayScale[i] = (uchar *) malloc(sizeof(uchar) * width);

    // Receive the submatrix from the proccess 0
    for (int i = 0; i < height; ++i)
      MPI_Recv(grayScale[i], width, MPI_UNSIGNED_CHAR, 0, 0,
         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Apply the filter
    applyFilter((void ***) &grayScale, fileType, filter, height, width);

    // Send the result back to proccess 0
    for (int i = 1; i < height - 1; ++i)
      MPI_Send(grayScale[i], width, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);

    // Free the memory of the matrix
    for (int i = 0; i < height; ++i)
      free(grayScale[i]);
    free(grayScale);

  } else if (!strcmp(fileType, "P6")) {
    grayScale = NULL;
    color = (pixel **) malloc(sizeof(pixel *) * height);
    for (int i = 0; i < height; ++i)
      color[i] = (pixel *) malloc(sizeof(pixel) * width);

    // Receive the submatrix from the proccess 0
    for (int i = 0; i < height; ++i)
      MPI_Recv(color[i], width, MPI_PIXEL, 0, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Apply the filter
    applyFilter((void ***) &color, fileType, filter, height, width);

    // Send the result back to proccess 0
    for (int i = 1; i < height - 1; ++i)
      MPI_Send(color[i], width, MPI_PIXEL, 0, 0, MPI_COMM_WORLD);

    // Free the memory of the matrix
    for (int i = 0; i < height; ++i)
      free(color[i]);
    free(color);
  }
}

int main(int argc, char *argv[]) {
  if (argc < 3)
    return -1;

  image input;
  int rank;
  int nProcesses;
  int height, width;
  char *fileType = (char *) malloc(sizeof(char) * 3);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

  // Create mpi data type for the structure of the pnm image
  const int structPixelItems = 3;
  int blocklengths[3] = {1, 1, 1};
  MPI_Datatype types[3] = {
    MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR
  };
  MPI_Datatype MPI_PIXEL;
  MPI_Aint offsets[3];

  offsets[0] = offsetof(pixel, red);
  offsets[1] = offsetof(pixel, green);
  offsets[2] = offsetof(pixel, blue);

  MPI_Type_create_struct(structPixelItems, blocklengths,
    offsets, types, &MPI_PIXEL);
  MPI_Type_commit(&MPI_PIXEL);

  // Read the image and send the info to the slave proccesses
  if (rank == 0) {
    readInput(argv[1], &input);
    int heightImage = atoi(input.height);
    int widthImage = atoi(input.width);

    height = heightImage / 3;
    height *= 3;
    height = height / nProcesses + 2;

    for (int i = 1; i < nProcesses; ++i) {
      MPI_Send(input.fileType, 3, MPI_CHAR, i, 0, MPI_COMM_WORLD);
      if (i != nProcesses - 1)
        MPI_Send(&height, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      else {
        int lastHeight = heightImage - (nProcesses - 1) * (height - 2);
        MPI_Send(&lastHeight, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      }
      MPI_Send(&widthImage, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }

  // Receive the needed data
  else {
    MPI_Recv(fileType, 3, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // For each filter call the operation functions
  if (rank == 0) {
    for (int i = 3; i < argc; ++i) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (!strcmp(argv[i], "bssembssem")) {
        int size = sizeof(bssembssem) / sizeof(bssembssem[0]);
        for (int j = 0; j < size; ++j) {
          MPI_Barrier(MPI_COMM_WORLD);
          operationProc0(&input, nProcesses, bssembssem[j], MPI_PIXEL);
        }
        continue;
      }
      operationProc0(&input, nProcesses, argv[i], MPI_PIXEL);
    }
    writeData(argv[2], &input);
  }

  // For each filter call the operation functions
  else {
    for (int i = 3; i < argc; ++i) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (!strcmp(argv[i], "bssembssem")) {
        int size = sizeof(bssembssem) / sizeof(bssembssem[0]);
        for (int j = 0; j < size; ++j) {
          MPI_Barrier(MPI_COMM_WORLD);
          operationSlaveProc(height, width, fileType, bssembssem[j], MPI_PIXEL);
        }
        continue;
      }
      operationSlaveProc(height, width, fileType, argv[i], MPI_PIXEL);
    }
  }

  free(fileType);
  MPI_Finalize();
  return 0;
}
