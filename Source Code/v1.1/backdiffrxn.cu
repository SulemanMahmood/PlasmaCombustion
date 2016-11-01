__global__ void reaction(double *a, int n)
{
  int tid = threadIdx.x;
  if (tid < 32)
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
}
