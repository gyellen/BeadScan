/* 
 *  This is a MEX file for MATLAB
 *  [xyt,frameIndex] = fastXYTFromList(photonList,dims,nFrames);  
 *  NOTE that dimension order for xyt is t,x,y
 *  dims should be (e.g.) [256 nLines nPixels]
 *
 *  Alternate call:
 *  fastXYTFromList(photonList,xyt,nFrames);
 *    USES existing xyt array (and its dimensions)
 */
 
#include "mex.h"
#include "matrix.h"
#include <Windows.h>
 
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
    byte *photonList;
    byte b;
    int dim_t, dim_lines, dim_pixels;
    int ix;
    size_t nBytes;
    UINT32_T *xyt_c, *xyt_src, *frameidx;                       /* pointer to dynamic data */
    mwSize ndim = 3;
    mwSize ndim_in;
    mwSize dims[3];
    mwSize elements;
    mwSize *xytdims;
    mwSize frameidx_dims[2] = {100,1};
    mxClassID arg2class;
    int k, iLine, iPixel, iFrame;
    int nFrames;
	
	if(nrhs != 3) {
		mexErrMsgIdAndTxt("MyToolbox:fastXYTFromList:nrhs",
						  "Two inputs required.");
	} 
    
    arg2class = mxGetClassID(prhs[1]);
    
    if (arg2class == mxUINT32_CLASS) {
        ndim_in = mxGetNumberOfDimensions(prhs[1]);
        xytdims = mxGetDimensions(prhs[1]);
        if(ndim_in != 3) {
            mexErrMsgIdAndTxt("MyToolbox:fastXYTFromList:argument2",
                    "Incoming XYT must be UINT32 and 3-dimensions (256,Lines,Pixels)");
        }
        // read the dimensions
        dim_t = (int) xytdims[0];
        dim_lines = (int) xytdims[1];
        dim_pixels = (int) xytdims[2];
        if(dim_t != 256) {
            mexErrMsgIdAndTxt("MyToolbox:fastXYTFromList:argument2",
                    "Incoming XYT must have dimension 1 = 256");
        }
        xyt_src = (UINT32_T *) mxGetData(prhs[1]);
    }
    else {
        // read the dimensions that are given directly by arg2
        dim_t      = (int) mxGetPr(prhs[1])[0];
        dim_lines  = (int) mxGetPr(prhs[1])[1];
        dim_pixels = (int) mxGetPr(prhs[1])[2];        
    }
    
	// create the output frameidx    
    /* Create a (nFrames+1)-by-1 mxArray */
	nFrames = (int) mxGetPr(prhs[2])[0];
	if (nFrames<1 || nFrames>100) {
            mexErrMsgIdAndTxt("MyToolbox:fastXYTFromList:argument3",
                    "Incoming nFrames must be between 1 and 100");		
	}
	frameidx_dims[0] = nFrames+1;
    plhs[1]  = mxCreateNumericArray(2, frameidx_dims, mxUINT32_CLASS, mxREAL);
	frameidx = (UINT32_T *) mxGetData(plhs[1]);

    // create the output array
    dims[0] = (mwSize) dim_t;
    dims[1] = (mwSize) dim_lines;
    dims[2] = (mwSize) dim_pixels;
    elements = dims[0] * dims[1] * dims[2];
    
    /* Create a local array with zeros */
    xyt_c = mxCalloc(elements, sizeof(UINT32_T));
    if (arg2class == mxUINT32_CLASS) {
        for(k=0; k<elements; k++) {
            // initialize the array with the starting values
            xyt_c[k] = xyt_src[k]; 
        }
    }
   
    /* get the input stream */
    photonList = (byte *)mxGetChars(prhs[0]);
    if(photonList == 0) {
        mexErrMsgIdAndTxt("MyToolbox:fastXYTFromList:photonList",
                "Byte array required.");
    }
    nBytes = mxGetNumberOfElements(prhs[0]);

    iLine = 0;
    iPixel = 0;
    iFrame = 0;
   
    for (k=0; k<nBytes; k++) {
        b = photonList[k];
        switch(b) {
            case 255 :
                iPixel = iPixel + 1;
                break;
            case 254 :
                iLine = iLine + 1;
                iPixel = 0;
                break;
            case 253 : 
                iFrame = iFrame + 1;
				if (iFrame <= nFrames) {
					frameidx[iFrame] = k+1;  // keep the frame index (at the i+1'th loc)
				}
                iPixel = 0;
                iLine = 0;
                break;
            default :
                ix = (iPixel * dim_lines + iLine) * dim_t + (int) b;
                if (ix<elements) {  // error check!
                    xyt_c[ix] = xyt_c[ix]+1;
                }
        }
    } 
   
    /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
    plhs[0] = mxCreateNumericArray(0, 0, mxUINT32_CLASS, mxREAL);

    /* Point mxArray to xyt_c */
    mxSetData(plhs[0], xyt_c);
    mxSetDimensions(plhs[0],dims,ndim);

    /* Do not call mxFree(dynamicData) because plhs[0] points to dynamicData */
    
    return;

}