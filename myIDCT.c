

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>



void printidct(int idct[640][480],int x,int y,FILE *fp)//function to save binary values to .pgm output file
{
	int i,j,temp,p,t;
	unsigned char c = 0x0;
	for(i=0;i<y;i++){//scan matrix from top to bottom
		for(j=0;j<x;j++){//scan each value of matrix from left to right
			temp=idct[j][i];//save value into temp
			c=0x0;//initialize c to 0x0
			for(p=0;p<8;p++){//convert temp from decimal to binary and save in c
				t=temp%2;
				if(t==1){
					c=c|1<<p;
				} else {
					c=c|0<<p;
				}
				temp=temp/2;
			}
			fwrite(&c,sizeof(c),1,fp);//write binary value to output.pgm file
		}
	}
}

void inversedct(double blockMatrix[8][8])//function to compute idct
{
	int u, v, x, y;
	double Cu, Cv, sum,idct,blockmatrixValue,tempdct[8][8];	
	for(y=0;y<8;y++)//scan block from top to bottom
	{
		for(x=0;x<8;x++)//scan each value from left to right
		{
			sum=0.0;
			for(v=0;v<8;v++){//compute value for Cu and Cv
				for(u=0;u<8;u++){
					if (u == 0) {
						Cu = 1.0 / sqrt(2.0);
					} else {
						Cu = 1.0;
					}	
					if (v == 0) {
						Cv = 1.0 / sqrt(2.0);
					} else {
						Cv = 1.0;
					} 
					//compute idct for each value of blockmatrix
					blockmatrixValue = blockMatrix[u][v];
				    idct = (double)(Cu*Cv*blockmatrixValue * cos(((2.0 * x)+1.0) * u * M_PI/16.0) * cos(((2.0 * y)+1.0) * v * M_PI/16.0));
					sum += idct;	
				}			
			}
			tempdct[x][y] = (0.25 * sum);//save idct value into tempdct matrix
			if(tempdct[x][y]>255) tempdct[x][y]=255;//crop value greater than 255 to 255
			if(tempdct[x][y]<0) tempdct[x][y]=0;//crop value less than 0 to 0
		}
	}
	for(u=0;u<8;u++)
	{
		for(v=0;v<8;v++){
			blockMatrix[v][u]=round(tempdct[v][u]);//store values from tempdct to blockMatrix
		}
	}
}

void iquant(double blockMatrix[8][8],int quantizationmatrix[8][8],float quantizationScaleValue)
{ //function to compute inverse of quantization
  int i,j;
  double mult;
	for(i=0;i<8;i++)//for each value in block
	{
		for(j=0;j<8;j++)
		{
			//subtract 127 and multiply by quantization coeff and qscale
			blockMatrix[j][i]-=127.0;
			mult=(double)blockMatrix[j][i]*((double)quantizationmatrix[j][i]*quantizationScaleValue);
			blockMatrix[j][i]=mult;
		}
	}
}

void izigzag(double block[8][8])
{//do inverse of zigzag reordering
	int i=0,j=0,c=0,r=0,nc=0,nr=0;
	double temp[8][8]={{0}};
	for(i=0;i<8;i++){//do inverse of zigzag reordering for top left corner of block
		if(i==0)//set temp[0][0]=block[0][0]
		{
			temp[r][c]=block[nr][nc];
			nr++;
		}else if((i%2)==0){//if i is even, read values sequentially from block and write out diagonal from top to bottom into temp
			r=0;c=i;//for matrix temp, set x co-ordinate to 0 and y co-ordinate to i
			while(c>=0)
			{
				if(nr==8) {//if y co-ordinate of block is 8, increase y co-ordinate by 1 and set x co-ordinate to 0
					nc++;
					nr=0;
				}
				temp[r][c]=block[nr][nc];//temp matrix is assigned corresponding block matrix value
				r++;//increase x co-ordinate by 1
				c--;//decrease y co-ordinate by 1
				nr++;//increase x co-ordinate of block matrix by 1
			}
		}
		else{//if i is odd, write out block values into diagonal (top to bottom) of temp matrix
			r=i,c=0;//for matrix temp,set x co-ordinate to i and y co-ordinate to 0
			while(r>=0)
			{
				if(nr==8){//if y co-ordinate of block is 8
					nc++;//increase y co-ordinate of block by 1
					nr=0;//set x co-ordinate of block to 0;
				}
				temp[r][c]=block[nr][nc];//temp matrix is assigned corresponding block matrix value
				c++;//increase y co-ordinate by 1
				r--;//decrease x co-ordinate by 1
				nr++;//increase x co-ordinate of block by 1
			}
		}
	}
	for(i=1;i<8;i++){//do inverse of zigzag reordering for the bottom right corner of block
		if((i%2)!=0)//if i is odd
		{
			r=i,c=7;//for matrix temp, set x co-ordinate to i and y co-ordinate to 7
			while(r<=7)
			{
				if(nr==8){//if y co-ordinate is 8
					nc++;//increase y co-ordinate of block by 1
					nr=0;//set x co-ordinate of block to 0
				}
				temp[r][c]=block[nr][nc];//temp is assigned corresponding block value
				c--;//decrease y co-ordinate by 1
				r++;//increase x co-ordinate by 1
				nr++;//increase x co-ordinate of block by 1
			}
		}
		else {//if i is even
			r=7;c=i;//for matrix temp,set x co-ordinate to 7 and y co-ordinate to i
			while(c<=7)
			{
				if(nr==8){//if y co-ordinate is 8
					nc++;//increase y co-ordinate of block by 1
					nr=0;//set x co-ordinate of block to 0
				}
				temp[r][c]=block[nr][nc];//temp is assigned corresponding block value
				r--;//decrease x co-ordinate by 1
				c++;//increase y co-ordinate by 1
				nr++;//increase x co-ordinate of block by 1
			}
		}
	}
	for(i=0;i<8;i++){
		for(j=0;j<8;j++){
			block[j][i]=temp[j][i];//save values from temp matrix to block matrix
		}
	}
}

int main(int argc, char *argv[])
{
	int quantization;
	int quantizationmatrix[8][8];
	float quantizationScaleValue;
	double block1[8][8],block2[8][8],block3[8][8],block4[8][8];//array to store gif matrix
	int idct[640][480],i,j,x,y,xunit,yunit,a,b,tempx,tempy;
	char xsize[10],ysize[10],symb[10],quantSV[10];
	FILE *fp1,*fp2,*fp3;
	if(argc!=4)// check if number of files specified is 3
	{
		printf("Error occured! Incorrect number of agruments are passed!");
		exit(1);
	}
	fp1=fopen(argv[1],"r");//open myDCT file (read only)
	fp2=fopen(argv[2],"r"); // opening the quantization file in read mode to access its elements
	fp3=fopen(argv[3],"wb");//open output .pgm file
	if(fp1==NULL)//check if myDCT file exists and can be opened
	{
		fprintf(stderr,"Can't open myDCT file");
		exit(1);
	}
	if(fp2  == NULL)  //check if quantization file exists
  	{
    	printf("Error opening %s for reading. Program terminated.", argv[2]);
    	return -1;
  	}
  	if(fp3  == NULL)  //check if output .pgm file can be opened
  	{
    	printf("Error opening %s for reading. Program terminated.", argv[2]);
    	return -1;
  	}
    for (i=0;i<8;i++) {//scan quantization file and save values into quantization matrix
        for(j=0;j<8;j++){
            fscanf(fp2, "%d", &quantization);
            quantizationmatrix[j][i] = quantization;
        }
    }
	for (i=0;i<4;i++)
	{
		fscanf(fp1,"%s",symb);
		if(i==1){ //get x size from myDCT file
			strcpy(xsize,symb);
			x=atoi(xsize);
			xunit=x/16;
		}
		if(i==2) {//get y size from myDCT file
			strcpy(ysize,symb);
			y=atoi(ysize);
			yunit=y/16;
		}
		if(i==3){//get qscale value from myDCT file
			strcpy(quantSV,symb);
			quantizationScaleValue=atof(quantSV);
		}
	}
	
	while(fscanf(fp1,"%d %d",&tempx,&tempy)!=EOF){//read each x y co-ordinate for a block
		for(i=tempy;i<(tempy+8);i++){//save block values into idct matrix
			for(j=tempx;j<(tempx+8);j++){
				fscanf(fp1,"%d",&a);
				idct[j][i]=a;
			}
		}
	}
	fprintf(fp3,"P5\n");//save magic marker into output.pgm file
	fprintf(fp3,"%d %d\n",x,y);//save x and y size into output.pgm file
	fprintf(fp3,"255\n");//save 255 to output.pgm file
	for(i=0;i<yunit;i++){//for each macroblock
		for(j=0;j<xunit;j++){
			for(a=0;a<8;a++){//divide each macroblock into 4 blocks
				for(b=0;b<8;b++){
					block1[b][a]=idct[j*16+b][i*16+a];//block 1 is the upper left block of macroblock
					block2[b][a]=idct[j*16+b+8][i*16+a];//block 2 is upper right block of macroblock
					block3[b][a]=idct[j*16+b][i*16+a+8];//block 3 is lower left block of macroblock
					block4[b][a]=idct[j*16+b+8][i*16+a+8];//block 4 is lower right block of macroblock
				}
			}
			izigzag(block1);//perform inverse zigzag for block 1
			izigzag(block2);//perform inverse zigzag for block 2
			izigzag(block3);//perform inverse zigzag for block 3
			izigzag(block4);//perform inverse zigzag for block 4
			iquant(block1,quantizationmatrix,quantizationScaleValue);//compute inverse quant values for block 1
			iquant(block2,quantizationmatrix,quantizationScaleValue);//compute inverse quant values for block 2
			iquant(block3,quantizationmatrix,quantizationScaleValue);//compute inverse quant values for block 3
			iquant(block4,quantizationmatrix,quantizationScaleValue);//compute inverse quant values for block 4
			inversedct(block1);//compute inverse of dct for block 1
			inversedct(block2);//compute inverse of dct for block 2
			inversedct(block3);//compute inverse of dct for block 3
			inversedct(block4);//compute inverse of dct for block 4
			for(a=0;a<8;a++){//for each of the 4 blocks
				for(b=0;b<8;b++){
					idct[j*16+b][i*16+a]=block1[b][a];//save block 1 values back to corresponding co-ordinates in idct
					idct[j*16+b+8][i*16+a]=block2[b][a];//save block 2 values back to corresponding co-ordinates in idct
					idct[j*16+b][i*16+a+8]=block3[b][a];//save block 3 values back to corresponding co-ordinates in idct
					idct[j*16+b+8][i*16+a+8]=block4[b][a];//save block 4 values back to corresponding co-ordinates in idct
				}
			}
		}
	}
	printidct(idct,x,y,fp3);//print out values from idct to output.pgm file
	fclose(fp1);//close mydct file
	fclose(fp2);//close quantization file
	fclose(fp3);//close output.pgm file
	return 0;
}


