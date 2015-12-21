
	
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>



void showblock(double block[8][8],FILE *fp)//function to save block coefficients to output file
{
	int i,j;
	for(i=0;i<8;i++){//for all values in block
		for(j=0;j<8;j++){
			fprintf(fp," %5d ",(int)block[j][i]);//save the 8*8 block to output file
		}
		fprintf(fp,"\n");
	}
}
void dct(double blockMatrix[8][8])//function to compute dct value for 8*8 block
{
	int u, v, i, j;
	double Cu, Cv, sum,dct,blockmatrixValue,tempdct[8][8];	
	for(u=0;u<8;u++)//for each u value, compute value of Cu
	{
		for(v=0;v<8;v++)//for each v value, compute value of Cv
		{
			if (u == 0) {//if u is 0 Cu is 1/sqrt(2) else it is 1
				Cu = 1.0 / sqrt(2.0);
			} else {
				Cu = 1.0;
			}	
			if (v == 0) {//if v is 0 Cv is 1/sqrt(2) else it is 1
				Cv = 1.0 / sqrt(2.0);
			} else {
				Cv = 1.0;
			} 
			sum = 0.0;	
			for (i = 0; i < 8; i++) {//scan matrix from top to bottom, left to right
				for (j = 0; j < 8; j++) {
					blockmatrixValue = blockMatrix[j][i];//save spacial value into blockmatrixValue and compute dct
			        dct = (double)(Cu*Cv*blockmatrixValue * cos(((2 * j)+1) * v * M_PI/16) * cos(((2 * i)+1) * u * M_PI/16));
					sum += (double)dct;//compute for entire block and save in sum		
				}				
			}
			tempdct[v][u] = 0.25 * sum;//the dct value for (v,u) is stored in tempdct matrix
		}
	}
	for(u=0;u<8;u++)
	{
		for(v=0;v<8;v++)
		{
			blockMatrix[v][u]=tempdct[v][u];//save values from tempdct matrix to original 8*8 block
		}
	}
}
void finalmatrix(double blockMatrix[8][8],int quantizationmatrix[8][8],float quantizationScaleValue)
{ //function for quantization of a 8*8 block
  int i,j;
  float division;
  for(i=0;i<8;i++)//scan from top to bottom, left to right
	{	
		for(j=0;j<8;j++)
		{
			division=blockMatrix[j][i]/(quantizationmatrix[j][i]*quantizationScaleValue);
			blockMatrix[j][i]=(int)round(division);//divide each value by quantization coefficient and qscale
			if(blockMatrix[j][i] > 128)//if value greater than 128, crop it to 128
			{
				blockMatrix[j][i] = 128;
			}
			else if(blockMatrix[j][i] < -127)//if value less than -127, crop it to -127
			{
				blockMatrix[j][i] = -127;
			}
			blockMatrix[j][i]+=127;
		}
	}
}

void zigzag(double block[8][8])
{//function to reorder 8*8 matrix in zigzag manner
	int i=0,j=0,c=0,r=0,nc=0,nr=0;
	double temp[8][8]={{0.0}};
	for(i=0;i<8;i++){//reorder upper left corner of matrix and save in 8*8 temp matrix
		if(i==0)
		{
			temp[nr][nc]=block[r][c];//block[0][0] is directly copied to temp[0][0]
			nr++;//increase x co-ordinate of temp matrix by 1
		}else if((i%2)==0){//when i is even, write out diagonal of block from bottom to top into temp
			r=0;c=i;// for block matrix set y co-ordinate to i and x co-ordinate to 0 
			while(c>=0)//while y is more than or equal to zero
			{
				if(nr==8) {//if x co-ordinate of temp matrix is 8, reset to 0 and increase y co-ordinate by 1
					nc++;
					nr=0;
				}
				temp[nr][nc]=block[r][c];//store block matrix value to temp matrix
				r++;//increase x co-ordinate by 1
				c--;//decrease y co-ordinate by 1
				nr++;//increase x co-ordinate of temp matrix by 1
			}
		}
		else{// when i is odd, write out diagonal from top to bottom
			r=i,c=0;//for block matrix, set x co-ordinate to i and y co-ordinate to 0
			while(r>=0)//while x is more than or equal to zero
			{
				if(nr==8){//if x co-ordinate of temp matrix is 8, reset to 0 and increase y co-ordinate by 1
					nc++;
					nr=0;
				}
				temp[nr][nc]=block[r][c];//store block matrix value to temp matrix
				c++;//increase y co-ordinate by 1
				r--;//increase x co-ordinate by 1
				nr++;//increase x co-ordinate of temp matrix by 1
			}
		}
	}
	for(i=1;i<8;i++){//perform zigzag reordering of bottom right corner of block matrix
		if((i%2)!=0)//for odd values of i, write out diagonal from top to bottom
		{
			r=i,c=7;//for block matrix, set x co-ordinate to i and y co-ordinate to 7
			while(r<=7)//while x co-ordinate is less than or equal to 7
			{
				if(nr==8){//if x co-ordinate of temp matrix is 8, reset to 0 and increase y co-ordinate by 1
					nc++;
					nr=0;
				}
				temp[nr][nc]=block[r][c];//store block matrix value to temp matrix
				c--;//decrease y co-ordinate by 1
				r++;//increase x co-ordinate by 1
				nr++;//increase x co-ordinate of temp matrix by 1
			}
		}
		else {//for even values of i, write out diagonal from bottom to top
			r=7;c=i;//for block matrix, set x co-ordinate to 7 and y co-ordinate to i
			while(c<=7)//while y co-ordinate is less than or equal to 7
			{
				if(nr==8){//if x co-ordinate of temp matrix is 8, reset to 0 and increase y co-ordinate by 1
					nc++;
					nr=0;
				}
				temp[nr][nc]=block[r][c];//store block matrix value to temp matrix
				r--;//decrease x co-ordinate by 1
				c++;//increase y co-ordinate by 1
				nr++;//increase x co-ordinate of temp matrix by 1
			}
		}
	}
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			block[j][i]=temp[j][i];//store value from temp matrix back to block matrix
		}
	}
}

int main(int argc, char *argv[])
{
	int gif[640][480];
	int quantization;
	int quantizationmatrix[8][8];
	float quantizationScaleValue;
	double block1[8][8],block2[8][8],block3[8][8],block4[8][8];//array to store gif matrix
	int i,j,x,y,xunit,yunit,a,b;
	char xsize[10],ysize[10],symb[10];
	FILE *fp1,*fp2,*fp3;
	if(argc!=5)// check if number of files specifies is 4
	{
		printf("Error occured! Incorrect number of agruments are passed!");
		exit(1);
	}
	fp1=fopen(argv[1],"r");//open pgm file (read only)
	fp2=fopen(argv[2],"r"); // open the quantization file in read mode to access its elements
	fp3=fopen(argv[4],"w");//open output file
	if(fp1==NULL)//check if pgm file exists and can be opened
	{
		fprintf(stderr,"Can't open pgm file");
		exit(1);
	}
	 if(fp2  == NULL)  // check if Quantization file can be opened
  	{
    	printf("Error opening %s for reading. Program terminated.", argv[2]);
    	return -1;
  	}
   for (i=0;i<8;i++) {//scan from top to bottom
        for(j=0;j<8;j++){//scan value of file from left to right
            fscanf(fp2, "%d", &quantization);//scan quantization value from file and save in quantizationmatrix;
            quantizationmatrix[j][i] = quantization;
        }
    }
	sscanf(argv[3], "%f", &quantizationScaleValue);//save quantization scale 
	for (i=0;i<4;i++)
	{
		fscanf(fp1,"%s",symb);//scan header of .pgm file
		if(i==1){ //scan x size of image
			strcpy(xsize,symb);
			x=atoi(xsize);
			xunit=x/16;
		}
		if(i==2) {//scan y size of image
			strcpy(ysize,symb);
			y=atoi(ysize);
			yunit=y/16;
		}
	}
	fgetc(fp1);//scan '\n' char
	for(i=0;i<y;i++){//scan file from top to bottom
		for(j=0;j<x;j++){//scan file from left to right
			gif[j][i]=fgetc(fp1);//save scanned file into gif matrix
		}
	}
	fprintf(fp3,"MYDCT\n");//save this to output file
	fprintf(fp3,"%d %d\n",x,y);//save x value and y value 
	fprintf(fp3,"%3f\n",quantizationScaleValue);//save quantization value to output file
	for(i=0;i<yunit;i++){//for each macroblock
		for(j=0;j<xunit;j++){
			for(a=0;a<8;a++){//for each block in macroblock
				for(b=0;b<8;b++){
					block1[b][a]=gif[j*16+b][i*16+a];//save coefficients from macroblock to block 1
					block2[b][a]=gif[j*16+b+8][i*16+a];//save coefficients from macroblock to block 2
					block3[b][a]=gif[j*16+b][i*16+a+8];//save coefficients from macroblock to block 3
					block4[b][a]=gif[j*16+b+8][i*16+a+8];//save coefficients from macroblock to block 4
				}
			}
			dct(block1);//compute dct for block 1 i.e. upper left 8*8 block of macroblock
			dct(block2);//compute dct for block 2 i.e. upper right 8*8 block of macroblock
			dct(block3);//compute dct for block 3 i.e. lower left 8*8 block of macroblock
			dct(block4);//compute dct for block 4 i.e. lower right 8*8 block of macroblock
			finalmatrix(block1,quantizationmatrix,quantizationScaleValue);//compute quantization for block 1
			finalmatrix(block2,quantizationmatrix,quantizationScaleValue);//compute quantization for block 2
			finalmatrix(block3,quantizationmatrix,quantizationScaleValue);//compute quantization for block3
			finalmatrix(block4,quantizationmatrix,quantizationScaleValue);//compute quantization for block 4
			zigzag(block1);// compute zigzag for block 1
			zigzag(block2);//compute zigzag for block 2
			zigzag(block3);//compute zigzag for block 3
			zigzag(block4);//computer zigzag for block 4
			fprintf(fp3,"%d %d\n",j*16,i*16);//write out x co-ordinate and y co-ordinate of top left coefficient of block 1
			showblock(block1,fp3);//save block 1 to output file
			fprintf(fp3,"%d %d\n",j*16+8,i*16);//write out x co-ordinate and y co-ordinate of top left coefficient of block 2
			showblock(block2,fp3);//save block 2 to output file
			fprintf(fp3,"%d %d\n",j*16,i*16+8);//write out x co-ordinate and y co-ordinate of top left coefficient of block 3
			showblock(block3,fp3);//save block 3 to output file
			fprintf(fp3,"%d %d\n",j*16+8,i*16+8);//write out x co-ordinate and y co-ordinate of top left coefficient of block 4
			showblock(block4,fp3);//save block 4 to output file
		}
	}
	fclose(fp1);//close pgm file
	fclose(fp2);//close quantization file
	fclose(fp3);//close output file
	return 0;
}


