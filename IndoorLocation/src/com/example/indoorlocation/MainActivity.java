package com.example.indoorlocation;

import java.io.BufferedOutputStream;
import com.example.indoorlocation.fftwork;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Timer;
import java.util.TimerTask;

import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.os.AsyncTask;
import android.os.Bundle;
import android.os.Environment;
import android.os.Handler;
import android.R.integer;
import android.R.string;
import android.app.Activity;
import android.content.Context;
import android.util.Log;
import android.view.Menu;
import android.view.View;
import android.view.View.OnClickListener;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

public class MainActivity extends Activity implements SensorEventListener
{

	private static final String TAG="RAWAUDIO";
	private boolean isRecording=false;
	private static int frequency=44100;
	String ip="192.168.1.100";// 服务端IP地址
	Handler handler=null;
	String msg;
	int c=0;
	int pc=0;
	int size=0;
	SensorManager sensorManager=null;
	float[] values;
	int leftorright=-1;

	@Override
	protected void onCreate(Bundle savedInstanceState)
	{
		sensorManager=(SensorManager)getSystemService(Context.SENSOR_SERVICE);
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_main);
		handler=new Handler();
		Button start=(Button)findViewById(R.id.start);
		Button stop=(Button)findViewById(R.id.stop);
		start.setOnClickListener(new OnClickListener() // 开始
		{
			public void onClick(View v)
			{
				isRecording=true;
				connectToServer cts=new connectToServer();
				cts.start();
			}
		});
		stop.setOnClickListener(new OnClickListener()
		{
			public void onClick(View v)
			{
				isRecording=false;
			}
		});
	}

	class connectToServer extends Thread
	{
		private String s;
		Socket sc;// 创建Socket连接
		DataOutputStream dout;
		int time;
		int end=0;

		public void run()
		{
			double[][] result=new double[2][2048];
			double[][] temp=new double[2][2048];
			int channelConfiguration=AudioFormat.CHANNEL_IN_STEREO;
			int audioEncoding=AudioFormat.ENCODING_PCM_16BIT;
			EditText x=(EditText)findViewById(R.id.edtx);
			EditText y=(EditText)findViewById(R.id.edty);
			EditText port=(EditText)findViewById(R.id.edtport);
			EditText edtip=(EditText)findViewById(R.id.edtip);
			ip=edtip.getText().toString();
			String sx=x.getText().toString();
			String sy=y.getText().toString();
			String sport=port.getText().toString();
			try
			{
				System.out.print("22");
				sc=new Socket(ip,Integer.parseInt(sport));
				dout=new DataOutputStream(sc.getOutputStream());

				int bufferSize=AudioRecord.getMinBufferSize(frequency,
						channelConfiguration,audioEncoding);
				short[] buffer=new short[bufferSize];
				AudioRecord audioRecord=new AudioRecord(
						MediaRecorder.AudioSource.MIC,frequency,
						channelConfiguration,audioEncoding,bufferSize);
				audioRecord.startRecording();
				double[] s1=new double[2050];
				double[] s2=new double[2050];
				int counter=0;
				int num=0;
				int end;
				while(isRecording)
				{
					int bufferReadResult=audioRecord.read(buffer,0,bufferSize);
					size=bufferReadResult;
					for(int i=0;i<4096;i++)
					{
						if(i%2==0)// 左声道
						{
							s1[i/2]=buffer[i];
							counter=i/2;
						}
						if(i%2==1)// 右声道
						{
							s2[counter]=buffer[i];
						}
					}
					fftwork.Accumulate(s1,s2,temp,result);
					temp=result;
					time++;
					if(time==10)
					{
						end=fftwork.Computing_TDOA(result[0],result[1]);
						if(end!=1999)
						{
							if(end>0)
							{
								dout.writeBytes("1"+" 方向:"+values[0]+""+" 坐标:"
										+"("+sx+","+sy+")\t\n");
								leftorright=1;
							}else
							{
								dout.writeBytes("0"+" 方向:"+values[0]+""+" 坐标:"
										+"("+sx+","+sy+")\t\n");
								leftorright=0;
							}
						}
						time=0;
						result=new double[2][2048];
						temp=new double[2][2048];
					}
				}
				audioRecord.stop();
				dout.close();
				sc.close();// 关闭Socket连接
			}catch(Throwable t)
			{
				
			}finally
			{
			}
		}

	}

	@Override
	protected void onResume()
	{
		super.onResume();
		sensorManager.registerListener(this,
				sensorManager.getDefaultSensor(Sensor.TYPE_ORIENTATION),
				SensorManager.SENSOR_DELAY_GAME);
	}

	Runnable runnableUI=new Runnable()
	{
		@Override
		public void run()
		{
			TextView txtleftorright=(TextView)findViewById(R.id.leftorright);
			txtleftorright.setText(""+leftorright);
		}
	};

	@Override
	public boolean onCreateOptionsMenu(Menu menu)
	{
		getMenuInflater().inflate(R.menu.main,menu);
		return true;
	}

	@Override
	public void onAccuracyChanged(Sensor arg0,int arg1)
	{
	}

	@Override
	public void onSensorChanged(SensorEvent event)
	{
		values=event.values;
		TextView txtdir=(TextView)findViewById(R.id.dir);
		txtdir.setText(values[0]+"");
	}
}
