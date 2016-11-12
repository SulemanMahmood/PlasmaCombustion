__constant__ k[29]

__global__ void reaction(double *species, int n)
{
  int tid = threadIdx.x;
  __shared__ s[928];
  s[tid] = species[tid];
  for (int t = 0; t < n; t++){
    if (tid < 32)
      s[tid] = -(k[0] + k[1] + k[2] + k[3] + k[4])*s[tid]*s[tid-896] + k[15]* + k[16]* + k[17]*
    else if (tid < 64)
    else if (tid < 96)
    else if (tid < 128)
    else if (tid < 160)
    else if (tid < 192)
    else if (tid < 224)
    else if (tid < 256)
    else if (tid < 288)
    else if (tid < 320)
    else if (tid < 352)
    else if (tid < 384)
    else if (tid < 416)
    else if (tid < 448)
    else if (tid < 480)
    else if (tid < 512)
    else if (tid < 544)
    else if (tid < 576)
    else if (tid < 608)
    else if (tid < 640)
    else if (tid < 672)
    else if (tid < 704)
    else if (tid < 736)
    else if (tid < 768)
    else if (tid < 800)
    else if (tid < 832)
    else if (tid < 864)
    else if (tid < 896)
    else if (tid < 928)
    __syncthreads();
  }
  species[tid] = s[tid];
}
