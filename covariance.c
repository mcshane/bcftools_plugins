/*  plugins/covariance.c -- XXX

    Copyright (C) 2015 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <math.h>
#include <getopt.h>

const char *about(void)
{
    return "Prints genotype genotype covariance matrix.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Print genotype covariance matrix\n"
        "Usage: bcftools +covariance [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Example:\n"
        "   bcftools +covariance in.vcf\n"
        "\n"
        "Definition:\n"
        "   C(x,y) = sum_{sites i} (x_i - 2*f_i)*(y_i - 2*f_i) /  (2 * f_i * (1 - f_i))\n"
        "\n";
}

bcf_hdr_t *in_hdr = NULL;
uint8_t *buf = NULL;
int nbuf = 0;   // NB: number of elements, not bytes
int ntags = 0, nsamples = 0;
double *cov, *dsg;

// C(x,y) = sum_{sites i} (x_i - 2*f_i)*(y_i - 2*f_i) /  (2 * f_i * (1 - f_i))
int site_covariance(bcf1_t *rec)
{
    int i, j, k, nret = bcf_get_genotypes(in_hdr,rec,(void**)&buf,&nbuf);
    if ( nret<0 ) return -1;
    double ac = 0, an = 0;

    nret /= rec->n_sample;
    int32_t *ptr = (int32_t*) buf;

    for (i=0; i<rec->n_sample; i++) dsg[i]=0;
    for (i=0; i<rec->n_sample; i++)
    {
        for (k=0; k<nret; k++)
        {
            if ( ptr[k]==bcf_int32_missing || ptr[k]==bcf_int32_vector_end || ptr[k]==bcf_gt_missing ) break;
            an++;
            if ( bcf_gt_allele(ptr[k]) ) { dsg[i]++; ac++; }
        }
        ptr += nret;
    }
    for (i=0; i<rec->n_sample; i++)
    {
        for (j=0; j<rec->n_sample; j++)
        {
            if (ac==0||ac==an) continue;
            cov[i*nsamples+j]+=(dsg[i]-2*ac/an)*(dsg[j]-2*ac/an)/(2*ac/an*(1-ac/an));
        }
    }
    return 0;
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int c;
    static struct option loptions[] =
    {
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?h",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }

    in_hdr = in;

    nsamples = bcf_hdr_nsamples(in_hdr);
    cov = (double*) calloc(nsamples*nsamples,sizeof(double));
    dsg = (double*) calloc(nsamples,sizeof(double));

    return 1;
}


bcf1_t *process(bcf1_t *rec)
{
    if (rec->n_allele!=2) return NULL; // biallelic sites only
    site_covariance(rec);
    return NULL;
}


void destroy(void)
{
    int i, j;
    printf("# COV, covariance\n");
    for (i=0; i<nsamples; i++)
        printf("\t%s", bcf_hdr_int2id(in_hdr,BCF_DT_SAMPLE,i));
    printf("\n");
    for (i=0; i<nsamples; i++)
    {
        printf("%s", bcf_hdr_int2id(in_hdr,BCF_DT_SAMPLE,i));
        for (j=0; j<nsamples; j++)
            printf("\t%.10f", cov[i<=j ? i*nsamples+j : j*nsamples+i]);
        printf("\n");
    }
    free(dsg);
    free(cov);
    free(buf);
}


