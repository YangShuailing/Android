package com.example.indoorlocation;

import java.util.Random;

import com.example.indoorlocation.Hamming;
import com.example.indoorlocation.filter_coff;

public class fftwork implements Anti_echo
{

	public static int MAX_FFT_SIZE=2048;
	public static int LF=8;
	public static int HF=372;
	public static int Max_D=62;
	public static int Max_H=30;
	public static int STESTNUMBER=16384;

	static void fft(double x[],double y[],int n,int sign)
	{
		int i,j,k,l,m=0,n1,n2;
		double c,c1,e,s,s1,t,tr,ti;
		for(j=i=1;i<16;i++)
		{
			m=i;
			j=2*j;
			if(j==n)
			{
				break;
			}
		}

		for(n1=n-1,j=0,i=0;i<n1;i++)
		{
			if(i<j)
			{
				tr=x[j];
				ti=y[j];
				x[j]=x[i];
				y[j]=y[i];
				x[i]=tr;
				y[i]=ti;
			}
			k=n/2;
			while(k<(j+1))
			{
				j=j-k;
				k=k/2;
			}
			j=j+k;
		}

		for(n1=l=1;l<=m;l++)
		{
			n1=2*n1;
			n2=n1/2;
			e=3.14159265359/n2;
			c=1.0;
			s=0.0;
			c1=Math.cos(e);
			s1=-sign*Math.sin(e);
			for(j=0;j<n2;j++)
			{

				for(i=j;i<n;i+=n1)
				{
					k=i+n2;
					tr=c*x[k]-s*y[k];
					ti=c*y[k]+s*x[k];
					x[k]=x[i]-tr;
					y[k]=y[i]-ti;
					x[i]=x[i]+tr;
					y[i]=y[i]+ti;
				}
				t=c;
				c=c*c1-s*s1;
				s=t*s1+s*c1;
			}
		}
		if(sign==-1)
		{
			for(i=0;i<n;i++)
			{
				x[i]/=n;
				y[i]/=n;
			}
		}
	}

	static int timedelay(double R[],int max_delay)

	{
		int k;
		double mean,mainpeak,hypopeak,tmp;
		int mainpeak_loc=0,hypopeak_loc;

		mean=0;
		mainpeak=-1;
		hypopeak=0;

		tmp=0;
		for(k=1;k<2*max_delay;k++) // 均值
		{
			tmp=tmp+R[k];
		}
		mean=tmp/(2.0*max_delay);

		for(k=1;k<2*max_delay;k++) /* 边界点不考虑 */
		{

			// if((R[k]>mainpeak)&&(R[k]>R[k-1])&&(R[k]>R[k+1]))
			if((R[k]>mainpeak)&&(R[k]>=R[k-1])&&(R[k]>=R[k+1]))
			{
				mainpeak=R[k];
				mainpeak_loc=k;
			}
		}

		if(mainpeak==-1)
			return -1000; /* 相关曲线中更本没有峰值 */

		hypopeak=mean; // 次峰赋初值为平均值
		hypopeak_loc=1;

		for(k=1;k<2*max_delay;k++)
		{
			if(k==mainpeak_loc)
				continue; // 不算主峰
			else if((hypopeak<R[k])&&(R[k]>R[k-1])&&(R[k]>R[k+1]))
			{
				hypopeak=R[k];
				hypopeak_loc=k;// 暂时尚没利用该信息
			}
		}

		if((hypopeak/mainpeak)<0.8)
			return mainpeak_loc-max_delay; /* 主峰足够尖锐，返回时延 */
		else
			return -2000; /* 由于次峰也很大，时延难以决定 */

	}

	static double[][] Accumulate(double s0[],double s1[],double temp[][],
			double result[][])
	{
		double[] real_f0=new double[2048];
		double[] image_f0=new double[2048];
		double[] real_f1=new double[2048];
		double[] image_f1=new double[2048];
		// double [] Re = new double[2048];
		// double [] Im = new double [2048];

		for(int k=0;k<MAX_FFT_SIZE;k++)
		{
			s0[k]=s0[k]*Hamming.Hamming[k];
			s1[k]=s1[k]*Hamming.Hamming[k];
		}

		for(int i=0;i<2048;i++)
		{
			real_f0[i]=s0[i];
			image_f0[i]=0;
		}

		fft(real_f0,image_f0,2048,1);

		for(int i=0;i<2048;i++)
		{
			real_f1[i]=s1[i];
			image_f1[i]=0;
		}
		fft(real_f1,image_f1,2048,1);

		for(int k=0;k<2048;k++)
		{
			result[0][k]=0;
			result[1][k]=0;
		}

		for(int i=LF;i<HF;i++)
		{
			result[0][i]=real_f0[i]*real_f1[i]+image_f0[i]*image_f1[i];
			result[1][i]=real_f0[i]*image_f1[i]*(-1)+image_f0[i]*image_f1[i];
		}

		for(int i=0;i<2048;i++)
		{
			result[0][i]=result[0][i]+temp[0][i];
			result[1][i]=result[1][i]+temp[1][i];
		}
		return result;
	}

	static int Computing_TDOA(double Re[],double Im[])
	{
		for(int i=LF;i<=HF;i++)
		{
			if(Re[i]!=0||Im[i]!=0)
			{
				Re[i]=Re[i]/Math.sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
				Im[i]=Im[i]/Math.sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
			}
		}

		// 实数的fft《---》共轭对称
		for(int k=LF;k<=HF;k++)
		{
			Re[2048-k]=Re[k];
			Im[2048-k]=-Im[k];
		}

		fft(Re,Im,2048,-1);

		double[] R01=new double[2048];
		for(int k=0;k<HF;k++)
		{
			R01[k]=0;
		}

		for(int k=-Max_D-1;k<0;k++)
		{
			R01[Max_D+k+1]=Re[2048+k];
		}
		for(int k=0;k<=Max_D+1;k++)
		{
			R01[k+Max_D+1]=Re[k];
		}

		int tao01=timedelay(R01,Max_D+1); // 垂直

		// System.out.println(-tao01-1);
		int ttdoa=-tao01-1;
		// if(ttdoa > 0 )
		// {
		// ttdoa =1;
		// }
		// else
		// {
		// ttdoa = 0;
		// }
		return ttdoa;
	}

	public static int TDOA(double s0[],double[] s1,double result)
	{
		double[] real_f0=new double[2048];
		double[] image_f0=new double[2048];
		double[] real_f1=new double[2048];
		double[] image_f1=new double[2048];
		double[] Re=new double[2048];
		double[] Im=new double[2048];

		for(int k=0;k<MAX_FFT_SIZE;k++)
		{
			s0[k]=s0[k]*Hamming.Hamming[k];
			s1[k]=s1[k]*Hamming.Hamming[k];
		}

		for(int i=0;i<2048;i++)
		{
			real_f0[i]=s0[i];
			image_f0[i]=0;
		}

		fft(real_f0,image_f0,2048,1);

		for(int i=0;i<2048;i++)
		{
			real_f1[i]=s1[i];
			image_f1[i]=0;
		}
		fft(real_f1,image_f1,2048,1);

		for(int k=0;k<2048;k++)
		{
			Re[k]=0;
			Im[k]=0;
		}

		for(int i=LF;i<HF;i++)
		{
			Re[i]=real_f0[i]*real_f1[i]+image_f0[i]*image_f1[i];
			Im[i]=real_f0[i]*image_f1[i]*(-1)+image_f0[i]*image_f1[i];
		}

		// //////////////////////////////////////////////////////////////////////
		// for (int i=LF;i<HF;i++)
		// {
		// Re[i] = Re[i] + real_f0[i] *real_f1[i] + image_f0[i]*image_f1[i];
		// Im[i] = Im[i] + real_f0[i] *image_f1[i]*(-1) + image_f0[i] *
		// image_f1[i];
		// }
		// }

		// //////////////////////////////////////////////////////////////////////
		// 谱白化
		for(int i=LF;i<=HF;i++)
		{
			if(Re[i]!=0||Im[i]!=0)
			{
				Re[i]=Re[i]/Math.sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
				Im[i]=Im[i]/Math.sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
			}
		}

		// 实数的fft《---》共轭对称
		for(int k=LF;k<=HF;k++)
		{
			Re[2048-k]=Re[k];
			Im[2048-k]=-Im[k];
		}

		fft(Re,Im,2048,-1);

		double[] R01=new double[2048];
		for(int k=0;k<HF;k++)
		{
			R01[k]=0;
		}

		for(int k=-Max_D-1;k<0;k++)
		{
			R01[Max_D+k+1]=Re[2048+k];
		}
		for(int k=0;k<=Max_D+1;k++)
		{
			R01[k+Max_D+1]=Re[k];
		}

		int tao01=timedelay(R01,Max_D+1); // 垂直

		// System.out.println(-tao01-1);
		int ttdoa=-tao01-1;
		if(ttdoa>0)
		{
			ttdoa=1;
		}else
		{
			ttdoa=0;
		}
		return ttdoa;
	}

	public static void main(String[] args)
	{
		int time=0;
		double[][] result=new double[2][2048];
		double[][] temp=new double[2][2048];
		for(int i=0;i<2;i++)
		{
			for(int j=0;j<2048;j++)
			{
				result[i][j]=0;
				temp[i][j]=0;
			}
		}
		double[] s0=new double[2048];
		double[] s1=new double[2048];
		int min=1;
		int max=100;
		Random ran=new Random();
		int shift=0;

		while(true)
		{
			for(int i=0;i<2048;i++)
			{
				s1[i]=max-ran.nextDouble()*(max-min);
			}

			for(int i=0;i<shift;i++)
			{
				s0[i]=0;
			}
			for(int i=shift;i<2048;i++)
			{
				s0[i]=s1[i-shift];
			}

			Accumulate(s0,s1,temp,result);
			temp=result;
			time++;
			if(time==100)
			{
				int end=Computing_TDOA(result[0],result[1]);
				System.out.println(end);

				time=0;
				result=new double[2][2048];
				temp=new double[2][2048];
		//		for(int i=0;i<2;i++)
		//		{
		//			for(int j=0;j<2048;j++)
		//			{
		//				result[i][j]=0;
		//				temp[i][j]=0;
		//			}
		//		}
				try
				{
					Thread.sleep(1000);
				}catch(InterruptedException e)
				{
					e.printStackTrace();
				}
			}
		}

		// computing_Tdoa=Accumulate(10);
		//
		// computing_Tdoa=Computed_TDOA(s0,s1);
		//
		//
		// computing_Tdoa=TDOA(s0,s1);

		// System.out.println(computing_Tdoa);

	}

}
