#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "efros_leung.h"

/**
 * @file    efros_leung_cmd.c
 * @brief   Main file for synthesis.
 * 
 * This file deals with I/O.
 * It defines the main function called from the command-line
 * with some arguments.
 */

/** @brief Struct of program parameters */
typedef struct
{
    int half_patch_sz;
    int out_sz;    
    int dimensions;
    float tolerance;
    char *input;
    char *output;
    char *synth_map_in;
    char *synth_map_out;

} program_argums;


static void usage(const char* name)
{
	printf("Efros Leung texture synthesis (accelerated version).\n");
    printf("Version 1.0 - 3 January 2013\n\n");
    printf("Usage: %s [options] <input file> <output file> <synthesis map sample file> <synthesis map output file>\n\n"
           "Only  PNG 8 bits images are supported.\n\n",name);
    printf("Options:\n");
    printf("  -p <number>  Half patch size");
    printf("  (default 4). Patches are of width 2p+1.\n");
	printf("  -t <number>  Tolerance parameter (default 0.1).\n");
    printf("  -s <number>  Output image size (default 128) \n");
	printf("  -d <number>  Number of PCA dimensions to keep, d <= (2p+1)p (default keep all).\n");
	printf("\n");
}

static void parse_arguments(program_argums *param, int argc, char *argv[])
{
    char *OptionString = NULL;
    char OptionChar;
    int i;


    if(argc < 1)
    {
        usage(argv[0]);
        exit(EXIT_SUCCESS);
    }


    /* loop to read parameters*/
    for(i = 1; i < argc;)
    {
        if(argv[i] && argv[i][0] == '-')
        {
            if((OptionChar = argv[i][1]) == 0)
            {
                error("Invalid parameter format.\n");
            }

            if(argv[i][2])
                OptionString = &argv[i][2];
            else if(++i < argc)
                OptionString = argv[i];
            else
            {
                error("Invalid parameter format.\n");
            }

            switch(OptionChar)
            {
                case 'p':
                    param->half_patch_sz = atoi(OptionString);
                    if(param->half_patch_sz < 1)
                    {
                        error("Invalid half patch size value.\n");
                    }
                    break;

                case 'd':
                    param->dimensions = atoi(OptionString);
                    if(param->dimensions < 1)
                    {
                        error("Invalid number of PCA dimensions.\n");
                    }
                    break;

                case 's':
                    param->out_sz = atoi(OptionString);                    
                    if(param->out_sz < 1)
                    {
                        error("Invalid output size.\n");
                    }
                    break;

                case 't':
                    param->tolerance = atof(OptionString);
                    if(param->tolerance < 0)
                    {
                        error("Invalid tolerance value.\n");
                    }
                    break;

                case '-':
                    usage(argv[0]);
                    exit(EXIT_FAILURE);

                default:
                    if(isprint(OptionChar))
                    {
                        fprintf(stderr, "Unknown option \"-%c\".\n",
                                OptionChar);
                        exit(EXIT_FAILURE);
                    } else
                        error("Unknown option.\n");
            }

        }
        else
        {
            if(!param->input)
                param->input = argv[i];
            else {
             if ( !param->output )
                param->output = argv[i];
             else { 
			  if ( !param->synth_map_in )    
                param->synth_map_in = argv[i];
              else   
                param->synth_map_out = argv[i];
			 }
			}
        }

        i++;
    }

    if( !param->input )
    {		
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }


    /* If parameters weren't set, set deafult parameters*/
    param->half_patch_sz = param->half_patch_sz>0 ? param->half_patch_sz : 4;
    param->out_sz = param->out_sz>0 ? param->out_sz : 128;    
    param->dimensions = param->dimensions>0 ? param->dimensions : (param->half_patch_sz*2+1)*param->half_patch_sz;    
    param->tolerance = param->tolerance>=0 ? param->tolerance : 0.1;
    param->output = param->output ? param->output : "output.png";
    param->synth_map_in = param->synth_map_in ? param->synth_map_in : "map.png";
    param->synth_map_out = param->synth_map_out ? param->synth_map_out : "copy_map.png";
    
}

/**
 * @brief   Main function.
 * 
 * The main function is used to:
 * @li get, check and initialize command-line parameters.
 * @li load and save the images.
 * @li launch the synthesis.
 */
int main(int argc, char ** argv)
{
  
  Cimage copy_map, map;  
  int init = 1;
	
  /*Initialize the structure param->* to -1 or null */
  program_argums param = {-1, -1, -1, -1, NULL, NULL, NULL, NULL};

  /*Parse command-line arguments*/
  parse_arguments(&param,argc,argv);
  
  Cimage w = load_png_image(param.input);  
  
  /* Verify if the input image is a true color image. Otherwise treat it as grayscale image. */
  /* This verification is only needed for the IPOL demo, since gray level input images are
   * converted to false color images */ 
  int is_grayScale = is_grayscale( w, w->ncol * w->nrow);  
  if ( is_grayScale == 1 )  {
    w->nchannels = 1;
  }
  
  /* Print algorithm information */
  print_algo_info(param.input, w->nrow, w->ncol, w->nchannels, param.half_patch_sz, param.tolerance, param.dimensions, param.out_sz);

  /* Launch the synthesis */
  Cimage v = efros_leung_synth( param.half_patch_sz, param.out_sz, param.tolerance, init, param.dimensions, w, &map, &copy_map);
  
  /* Write the output images */
  write_png_image(v, param.output); 
  write_png_image(map, param.synth_map_in);
  write_png_image(copy_map, param.synth_map_out);  
 
  /* Exit */
  delete_Cimage(w);
  delete_Cimage(v);
  return EXIT_SUCCESS;

}
