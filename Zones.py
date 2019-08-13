#这个版本利用张韵豪推到的数学公式，修订了Liu Qing等人在β角计算上的错误
global Log
Log=[]

from sys import exit
from shutil import copyfile
import os
import time

TIME=str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

with open("comment.txt","a") as f:
	f.write(str(TIME)+'\n')
	f.write('Program started\n\n')

TIME=TIME.replace(':','_')

File_Path = os.getcwd()+str('\saved_files\\')+TIME
#os.makedirs(File_Path)	#创建文件夹

def goodbye():
	with open("comment.txt","a") as f:
		f.write('LOG: '+str(Log)+'\n')
		f.write('*********************\n')
	sys.exit(0)

def leave_comment():	#弹出feedback界面
	def vset1():		#当yes的选项被点选时的函数
		global answer
		answer=str('right')

	def vset2():		#当NO被点选时
		global answer
		answer=str('wrong')
   
	def save_comment():	 #点击submit时，保存Y/N和评论
		def destroyall():
			cmt2.destroy()
			cmt.destroy()

		comment=comment_entry.get("1.0",'end-1c')   #这里加入存档指令	 
		with open("comment.txt","a") as f:
			f.write(str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))+'\nComment: '+answer+'___(')			 #格式：right/wrong_abcd
			f.write(comment)
			f.write(')\n\n')


		cmt2=tkinter.Tk()
		cmt2.title('Thank you')		 #存档结束后弹窗
		cmt2.geometry('500x150')
		tkinter.Label(cmt2,text='\n').pack()
		tkinter.Label(cmt2,text='Your words have been recorded. Thank you for helping us!').pack()
		tkinter.Label(cmt2,text='\n').pack()
		Button(cmt2,text='Done!',command=destroyall).pack() #关闭
		Log.append('X1')
		cmt2.mainloop()
	
	cmt=tkinter.Tk()		#feedback的主界面
	cmt.title('comment')
	cmt.geometry('700x540')

	tkinter.Label(cmt,text='\nDid it give a right answer?').pack()  #第一个单选题
	v=StringVar()
	rb1=Radiobutton(cmt, text='Yes, the answer was correct.',variable=v,value='1',command=vset1)
	rb2=Radiobutton(cmt, text='No, the answer was incorrect.',variable=v,value='2',command=vset2)
	rb1.pack()
	rb2.pack()
	rb2.select()		#保证两个选项是同时点选的
	rb2.deselect()

	
	tkinter.Label(cmt,text='Please comment in the entry below\n').pack()	
	comment_entry=Text(cmt)	 #Text widget用于输入框
	comment_entry.pack()
	tkinter.Label(cmt,text='\n').pack()
	Button(cmt,text='Submit!',command=save_comment).pack()
	Log.append('X')
	cmt.mainloop()

while True:			#保证循环，直到用户自主关闭

	def distance(a,b):	#已知两点坐标，算距离
		return ((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5

	try:
	#尝试，有错报错
		global state
		state=1
		#state状态参数：1代表定义函数阶段，加上起始输入hkl录入。这部分成功之后，会直接进入2状态而不回头
		#			   2代表图上选点阶段，通过读取步骤1内输入的图片，用户输入点

		while state==1:
			#在这一版本中，通过点击三个点来标定
			import math
			import time
			import numpy as np
			from scipy import stats
			from PIL import Image
			from pylab import *

			#比例尺sacle=?A/pix
			#scale=385

			###########################################
			import tkinter
			from tkinter import Button
			from tkinter import Label
			from tkinter import StringVar
			from tkinter import Entry
			from tkinter import Pack
			from tkinter import Text
			from tkinter import Radiobutton
			from tkinter import IntVar
			from tkinter import StringVar

			import tkinter.filedialog

			#################################################
			##########################################
			#输入相机常数
			def xiangji():
				global pixel_length
				pixel_length=e.get()
				pixel_length=float(pixel_length)
				win_xiangji.destroy()

			global win_xiangji
			global aconst
			win_xiangji=tkinter.Tk()
			win_xiangji.title('Pixel length')
			win_xiangji.geometry('300x180')
			aconst=StringVar()
			aconst.set('0.0001536')		   #缺省值
			tkinter.Label(win_xiangji,text='\nPlease check pixel size, (1/A)').pack()
			e=Entry(win_xiangji,textvariable=aconst,width=10)
			e.pack()

			global wave_length
			wave_length=0.02508				#输入波长，对于电压200V的TEM而言，波长固定为0.02508埃

			tkinter.Label(win_xiangji,text='\n').pack()
			Button(win_xiangji,text='done!',command=xiangji).pack()
			Button(win_xiangji,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)
			win_xiangji.mainloop()

			#**************************************#
			########################################

			##########################################
			#选取标定的电子衍射图片,然后Laue或者hkl文件###
			#*****************************************#

			def win_choose_image():
				global filename1
				global im
				filename1=tkinter.filedialog.askopenfilename()
				im=Image.open(filename1)
				button_hkl.pack()
				tkinter.Label(root,text='').pack()
				button_laue.pack()
		
			def choose_hkl():
				global filename2,ftext,A,B,C,X,Y,Z

				filename2=tkinter.filedialog.askopenfilename()
				############################################
				##read the input hkl file

				f=open(filename2)
				ftext=f.readlines()
				
				for i in range(2):
					ftext.pop(0)

				##########################
				#get ABC XYZ
				
				ftem=[]
				f=ftext[0].split(' ')
				for i in range(len(f)):
					if f[i]=='' or f[i]==';' or f[i]=='=' or f[i]=='Cell':
						pass

					else :
						ftem.append(f[i].replace('\n',''))

				A=float(ftem[0])
				B=float(ftem[1])
				C=float(ftem[2])
				X=float(ftem[3])*math.pi/180
				Y=float(ftem[4])*math.pi/180
				Z=float(ftem[5])*math.pi/180

				ftext.pop(0)
				ftext.pop(0)

				##########################
				#get hkl & d
				##########################

				#检查输入没有问题
				global hkl_ref,d_ref
				hkl_ref=[]
				d_ref=[]

				hkl_ref.append([1,0,0])	 #debug:HKL文件中不会包含基矢量
				d_ref.append(A)
				hkl_ref.append([0,1,0])
				d_ref.append(B)
				hkl_ref.append([0,0,1])
				d_ref.append(C)
				
				for i in range(len(ftext)):
					f=ftext[i].split(' ')
					hkl_ref.append([])
					k=0
					for j in range(len(f)):
						if f[j]=='' :
							pass
						elif k<3:
							hkl_ref[i+3].append(int(f[j]))
							k+=1
						elif k==3:
							d_ref.append(float(f[j]))
							k+=1

				#######################

				tkinter.Label(root,text='\n').pack()
				button_zone_zone.pack()

			#############################################################
			def click_center():
				l=im.size[0]	#读图片取长，图片为im
				w=im.size[1]	#读图片取宽
				s=int(min(l,w)) #两者取小
				global S
				S=int(0.3*s)	#设置局部放大倍率

				def strange():
					win_click_center.destroy()

					mng = plt.get_current_fig_manager()	 #让其全屏显示
					mng.window.state('zoomed')

					#展示图片，未放大，然后在图片上点一下#######
					#***************************************#
					imshow(im)
					zoom_point=ginput(1)	#一下
					close()
					#***************************************#
					#图片关闭################################
					
					zoom_center=zoom_point[0]   #先以第一次点的粗略位置为中心，放大确定衍射中心
					
					box_pre=(zoom_center[0]-S,zoom_center[1]-S,zoom_center[0]+S,zoom_center[1]+S)
					imzoomed_pre=im.crop(box_pre)   #框定范围，此时imzoomed_pre为目标图片，原图片为im

					mng=plt.get_current_fig_manager()	 #让其全屏显示
					mng.window.state('zoomed')

					#展示已放大的图片，然后在图片上点两下确定衍射中心####
					#************************************#
					imshow(imzoomed_pre)				
					sym_point_pre=ginput(2) #点两下
					close()
					#************************************#
					#图片关闭##############################		 

					sym_point=[]	   #将局部坐标转化全局坐标
					for i in sym_point_pre:
						sym_point.append([i[0]-S+zoom_center[0],i[1]-S+zoom_center[1]])

					global center 
					center=[0.5*(sym_point[0][0]+sym_point[1][0]),0.5*(sym_point[0][1]+sym_point[1][1])]

					global box
					box=(center[0]-S,center[1]-S,center[0]+S,center[1]+S)  #以衍射中心为中心，设置局部放大区域
					
					global xx,yy,point_A,point_B,OA_plot,OB_plot
					xx=50*(2)**0.5*math.cos(40.59/180*math.pi)	 
					yy=50*(2)**0.5*math.sin(40.59/180*math.pi)
					point_A=[center[0]+xx,center[1]-yy]
					point_B=[center[0]+xx,center[1]+yy]			 #OA是x,OB是y
					OA_plot=[xx,-yy]
					OB_plot=[xx,yy] ##########################


				############################
				#**************************#
				win_click_center=tkinter.Tk()
				win_click_center.title('Second Step')
				win_click_center.geometry('400x150')
				tkinter.Label(win_click_center,text='\nPlease click on the image to zoom in, then select the center. ').pack()
				tkinter.Button(win_click_center,text='OK',command=strange).pack()
				Button(win_click_center,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)

				Log.append('c')
				win_click_center.mainloop()	 #弹框结束
				#**************************#
				############################

			#############################################################

			def laue():
				win_holder.destroy()	#getinputdata获得样品杆角度，触发laue，于是结束上层窗口

				click_center()			#输入中心是必要的

				imshow(im)
				#plot(S,S,'r*') #可选的其他形状 o . s + x
				circle_points=ginput(40)		 
				close()

				##计算拟合圆的数据
				circle_x = []
				circle_y = []#这两个集合用来储存散点的x、y坐标

				for i in circle_points:
					circle_x.append(i[0])
					circle_y.append(i[1])
					
				#接下来用最小二乘法拟合得到圆
				laue_N = len(circle_x)
				sum_laue_x=0
				sum_laue_y=0
				sum_laue_xy=0
				sum_laue_x2=0
				sum_laue_y2=0
				sum_laue_x3=0
				sum_laue_y3=0
				sum_laue_x2y=0
				sum_laue_xy2=0
				i=0
				while i < len(circle_x):
					sum_laue_x += circle_x[i]
					sum_laue_y += circle_y[i]
					sum_laue_xy += circle_x[i]*circle_y[i]
					sum_laue_x2 += circle_x[i]**2
					sum_laue_y2 += circle_y[i]**2
					sum_laue_x3 += circle_x[i]**3
					sum_laue_y3 += circle_y[i]**3
					sum_laue_x2y += (circle_x[i]**2)*circle_y[i]
					sum_laue_xy2 += circle_x[i]*(circle_y[i]**2)
					i += 1

				laue_C = laue_N*sum_laue_x2-sum_laue_x*sum_laue_x
				laue_D = laue_N*sum_laue_xy-sum_laue_x*sum_laue_y
				laue_E = laue_N*(sum_laue_x3+sum_laue_xy2)-sum_laue_x*(sum_laue_x2+sum_laue_y2)
				laue_G = laue_N*sum_laue_y2-sum_laue_y*sum_laue_y
				laue_H = laue_N*(sum_laue_y3+sum_laue_x2y)-sum_laue_y*(sum_laue_x2+sum_laue_y2)

				laue_a = (laue_H*laue_D-laue_E*laue_G)/(laue_C*laue_G-laue_D**2)
				laue_b = (laue_H*laue_C-laue_E*laue_D)/(laue_D**2-laue_G*laue_C)
				laue_c = -(sum_laue_x2+sum_laue_y2+laue_a*sum_laue_x+laue_b*sum_laue_y)/laue_N

				#得到Laue圆的信息，圆心在(laue_A,laue_B)，半径为laue_R
				laue_A = -0.5*laue_a
				laue_B = -0.5*laue_b
				laue_R = 0.5*(laue_a**2+laue_b**2-4*laue_c)**0.5
				###################################################################
				#一些必要的输入，在这里提前放入

				#在这里统一所有变量，直接与张韵豪的论文公式对应
				point_O=center
				point_Ol=[laue_A,laue_B]

				#计算公式需要的角度
				len_OA=distance(point_O,point_A)
				len_OY=distance(point_O,point_B)
				len_OOl=distance(point_O,point_Ol)
				len_OlA=distance(point_Ol,point_A)
				len_OlY=distance(point_Ol,point_B)
				cos_YOOl=(len_OY**2+len_OOl**2-len_OlY**2)/(2*len_OY*len_OOl)
				cos_AOOl=(len_OA**2+len_OOl**2-len_OlA**2)/(2*len_OA*len_OOl)

				#计算转角公式需要的量
				laue_x=laue_R*cos_AOOl
				laue_y=laue_R*cos_YOOl
				laue_z=-1/(wave_length*pixel_length)
				laue_w=-(laue_z**2+laue_R**2)**(0.5)

				#计算laue需要的转角
				laue_OA_rotate=arcsin((laue_z*math.sin(x1)+laue_y*math.cos(x1))/laue_w)-x1
				laue_OB_rotate=-arctan(laue_x/(math.sin(x1)*laue_y-math.cos(x1)*laue_z))

				laue_OA_result=(laue_OA_rotate+x1)/math.pi*180
				laue_OB_result=(laue_OB_rotate+y1)/math.pi*180

				global state
				state=1

				win_show_laue=tkinter.Tk()
				win_show_laue.title('Laue Circle Result')
				win_show_laue.geometry('400x180')
				tkinter.Label(win_show_laue,text='\nTo tilt the crystal to a proper axis, tilt X and Y to: \n').pack()
				tkinter.Label(win_show_laue,text='Final X rod position: '+str(round(laue_OA_result,2))).pack()
				tkinter.Label(win_show_laue,text='Final Y rod position: '+str(round(laue_OB_result,2))+'\n').pack()
				tkinter.Button(win_show_laue,text='OK',command=win_show_laue.destroy).pack()
				win_show_laue.mainloop()


			def lauewindow():
				global win_holder
				root.destroy()
				win_holder=tkinter.Tk()
				win_holder.title('Special for Laue Circle')
				win_holder.geometry('400x200')

				tkinter.Label(win_holder,text='\nPlease give the current tilt angles of x & y (0-180)\n\n\n').pack()
				etiltx=Entry(win_holder,width=5)
				etiltx.place(x=150,y=40)
				etilty=Entry(win_holder,width=5)
				etilty.place(x=200,y=40)

				tkinter.Label(win_holder,text='Please click on the center,then 20 independent spots.').pack()

				def get_inputdata():
					global x1,y1
					x1=float(etiltx.get())*math.pi/180
					y1=float(etilty.get())*math.pi/180
					laue()
				
				button=tkinter.Button(win_holder,text='Done!',command=get_inputdata)
				button.pack()
				Button(win_holder,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)

				win_holder.mainloop()

			root=tkinter.Tk()
			root.title('First Step')
			root.geometry('600x420')
			tkinter.Label(root,text='\nPlease give the input files: (image) , (.txt) ').pack()

			button_choose_img=tkinter.Button(root,text='Choose the image',command=win_choose_image)
			button_choose_img.pack()

			button_hkl=tkinter.Button(root,text='Choose the hkl file for zone to zone calculation',command=choose_hkl)
			tkinter.Label(root,text='\n').pack()

			button_laue=tkinter.Button(root,text='	Or start Laue Circle refinement	',command=lauewindow)

			def close_root_click_center():
				global state
				state=2
				root.destroy()
				click_center()
				

			button_zone_zone=tkinter.Button(root,text='Calculate zone to zone angle',command=close_root_click_center)

			Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)
			tkinter.Label(root,text='').pack(side=tkinter.BOTTOM)
			Button(root,text='Any inconvience? Submit feedback',command=leave_comment).pack(side=tkinter.BOTTOM)

			Log.append('b')
			root.mainloop()

			#选取图片和hkl的UI窗口关闭**********************************#
			#############################################
			###########################################	
			
			#第一波输入已完成，创建一些函数
			#############################################
			   
			#funcition COMPARE:a is an item and b is a list, find the item closest to a and return its index in b
			def compare(a,b):
				x=[]
				for i in a:
					y=[]
					for j in range(len(b)):
						y.append(abs(b[j]-i))
					c=y.index(min(y))
					x.append(c)
				return x

			def find_close(a,b):		#a为数值，b为列表，返回
				diff=[]
				x1=[]
				x2=[]
				k1=0
				k2=0
				k3=0
				for i in b:
					diff.append([])
					# k=0
					for j in i:
						diff[k1].append(abs(a-j))
					k1=k1+1
				while k2<len(diff):
					   # print(k2)
						#if# diff[k2]==min(diff):
					while k3<len(diff[k2]):
						if diff[k2][k3]==min(diff[k2]):
							x2.append(k3)
							x1.append(k2)
							break
						k3=k3+1
					k3=0
					k2=k2+1
						   
				return [x1,x2]

			##To calculate angle between three spots
			def measureangle(center,AA,BB):
				a=[]
				b=[]
				a.append(AA[0]-center[0])
				a.append(AA[1]-center[1])
				b.append(BB[0]-center[0])
				b.append(BB[1]-center[1])
				x1=math.sqrt(a[0]**2+a[1]**2)
				x2=math.sqrt(b[0]**2+b[1]**2)
				x3=math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
				y=(x1**2+x2**2-x3**2)/(2*x1*x2)
				ang=round(math.acos(y),5)
				return ang  

			#求晶面夹角,input:[a,b,c],[d,e,f],d1,d2
			def plane_angle(hkl1,hkl2,d1,d2):
				h1=hkl1[0]
				k1=hkl1[1]
				l1=hkl1[2]
				h2=hkl2[0]
				k2=hkl2[1]
				l2=hkl2[2]
				s11=(B*C*math.sin(X))**2
				s22=(A*C*math.sin(Y))**2
				s33=(A*B*math.sin(Z))**2
				s12=A*B*C*C*(math.cos(X)*math.cos(Y)-math.cos(Z))
				s23=A*A*B*C*(math.cos(Y)*math.cos(Z)-math.cos(X))
				s13=A*B*B*C*(math.cos(Z)*math.cos(X)-math.cos(Y))
				V=A*B*C*(1-math.cos(X)**2-math.cos(Y)**2-math.cos(Z)**2+2*math.cos(X)*math.cos(Y)*math.cos(Z))**0.5
				re=math.acos(d1*d2/(V**2)*(s11*h1*h2+s22*k1*k2+s33*l1*l2+s23*(k1*l2+k2*l1)+s13*(l1*h2+l2*h1)+s12*(h1*k2+h2*k1)))
				return round(re,5)

			def ortho_angle(a,b):#vector a and b are in the orthonormal axis
				m=(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/((a[0]**2+a[1]**2+a[2]**2)*(b[0]**2+b[1]**2+b[2]**2))**0.5
				return math.acos(abs(m))

			def trans_zone(a):#将任意晶系的晶向矢量投射进立方晶系
				x1=1
				y1=0
				z1=0
				x2=math.cos(Z)
				y2=math.sin(Z)
				z2=0
				x3=math.cos(Y)
				y3=abs(math.cos(X))+abs(math.cos(Y)/math.cos(Z))/math.sin(Z)-abs(math.cos(Y)*tan(Z))
				z3=(abs(1-x3**2-y3**2))**0.5
				p=[a[0]*x1+a[1]*x2+a[2]*x3,a[0]*y1+a[1]*y2+a[2]*y3,a[0]*z1+a[1]*z2+a[2]*z3]
				p[0]=p[0]*A
				p[1]=p[1]*B
				p[2]=p[2]*C
				return p

			def zone_angle(HKL1,HKL2):
				m1=trans_zone(HKL1)
				m2=trans_zone(HKL2)
				return ortho_angle(m1,m2)

			def delpnzone(a):   #去除正负重复，如[011]和[0-1-1]，再把[0-1-1]变为[011]
				b=[]
				c=[]
				for i in a:
					if [i[0],i[1],i[2]] in b or [-i[0],-i[1],-i[2]] in b:
						pass
					else:
						b.append(i)
				for i in b:
					if (i[0]==abs(i[0]) and i[1]==abs(i[1])) or (i[0]==abs(i[0]) and i[2]==abs(i[2])) or (i[1]==abs(i[1]) and i[2]==abs(i[2])):
						c.append(i)
					else:
						c.append([-i[0],-i[1],-i[2]])
				return c

			def cal_positive(list):
					j=0
					for i in list:
							if i>0:
									j=j+1
							if i<0:
									j=j-1
							if i==0:
									j=j+0.5
					return j
						
			#化简[a,b,c],如[4,0,0]变成[1,0,0]
			def simplify(list):
				list1=[]
				if cal_positive(list)<=0:
						for i in list:
								list1.append(-i)
				else:
						for i in list:
								list1.append(i)

				sm=max(list1)
				re=[]
				re.append(list1[0])
				re.append(list1[1])
				re.append(list1[2])
				while sm>0:
					if re[0]%sm==0 and re[1]%sm==0 and re[2]%sm==0:
						re[0]=int(re[0]/sm)
						re[1]=int(re[1]/sm)
						re[2]=int(re[2]/sm)
					sm=sm-1
				return re

			#求与[a,b,c],[d,e,f]垂直的向量
			def ortho(list):
				a1=list[0]
				a2=list[1]
				m1=[a1[1]*a2[2]-a1[2]*a2[1],a1[2]*a2[0]-a1[0]*a2[2],a1[0]*a2[1]-a1[1]*a2[0]]
				return m1

			def delrepeat(list):
				b=[]
				for i in range(len(list)):
					if list[i] in b:
						pass
					else:
						b.append(list[i])
				return b

			def hklfamily(a):
				b=[]
				b.append([a[0],a[1],a[2]])
				b.append([(-1)*a[0],a[1],a[2]])
				b.append([a[0],(-1)*a[1],a[2]])
				b.append([a[0],a[1],(-1)*a[2]])
				b.append([(-1)*a[0],(-1)*a[1],a[2]])
				b.append([a[0]*(-1),a[1],(-1)*a[2]])
				b.append([a[0],(-1)*a[1],(-1)*a[2]])
				b.append([a[0]*(-1),a[1]*(-1),(-1)*a[2]])

				b=delrepeat(b)
				return b

			def hkilfamily(a):	  ##########六方or三方
				b=[]
				b.append([a[0],a[1],a[2],a[3]])
				b.append([a[2],a[0],a[1],a[3]])
				b.append([a[1],a[2],a[0],a[3]])
				b.append([a[1],a[0],a[2],(-1)*a[3]])
				b.append([a[0],a[2],a[1],(-1)*a[3]])
				b.append([a[2],a[1],a[0],(-1)*a[3]])
				b.append([(-1)*a[0],(-1)*a[1],(-1)*a[2],(-1)*a[3]])
				b.append([(-1)*a[2],(-1)*a[0],(-1)*a[1],(-1)*a[3]])
				b.append([(-1)*a[1],(-1)*a[2],(-1)*a[0],(-1)*a[3]])
				b.append([(-1)*a[1],(-1)*a[0],(-1)*a[2],a[3]])
				b.append([(-1)*a[0],(-1)*a[2],(-1)*a[1],a[3]])
				b.append([(-1)*a[2],(-1)*a[1],(-1)*a[0],a[3]])

				b=delrepeat(b)
				return b

			#################################  

			#####################################################
			while state==2:
				#必要的时候留存文件，但是调试的时候没有必要
				#copyfile(filename1,File_Path+'\\pattern')   #若不报错，则复制
				#copyfile(filename2,File_Path+'\\hkl')

				##################################
				#********************************#
				#弹窗，请点五对点
				root=tkinter.Tk()
				root.title('Third Step')
				root.geometry('400x150')
				tkinter.Label(root,text='\nPlease click on 5 pairs of independent spots.').pack()
				tkinter.Button(root,text='OK',command=root.destroy).pack()
				Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)

				Log.append('d')
				root.mainloop()	 #弹框结束
				#********************************#
				##################################

				mng = plt.get_current_fig_manager()
				mng.window.state('zoomed')

				#弹出放大后的图片，点十对点##########
				#********************************#
				imzoomed=im.crop(box)
				imshow(imzoomed)
				plot(S,S,'r*') #可选的其他形状 o . s + x
				ten_points_pre=ginput(10)		 #点10下，两对两对点
				close()
				#********************************#
				#图片关闭##########################
				
				#############################
				ten_points=[]
				for i in ten_points_pre:
					ten_points.append([i[0]-S+center[0],i[1]-S+center[1]])	  #整理，获得十个点的图上坐标在tenpoints
				#############################

				fin_yrw=[]	  #五组标定结果放进yrw
				showangle=[]	#展示每组的角度
				showhkl=[]
				showd_cal=[]	#展示每组的d
				showd_real=[]

				def index_zone(three_points,multi1,multi2):	#将选点比对d值和夹角进行标定，然后叉乘两晶面输出带轴的函数
					global fin_zone1,HKL1,HKL2,HKL3,pre_direction1,pre_direction2		
					global angle_plot1		#测量倒异空间中的夹角
					angle_plot1=measureangle(center,three_points[0],three_points[1])

					r=[]
					r.append(distance(three_points[0],center)*multi1)
					r.append(distance(three_points[1],center)*multi2)
					r.append(distance(three_points[2],center)*multi1)	#r代表图中的像素距离

					#从像素距离r转换成晶面间距d_cal，real distance = length on the map * pixel_length
					d_cal=[]	#测量值
					for i in r:
						d_cal.append(1/(i*pixel_length))

					#比较点中的晶面距离和参考中的数据
					p1=compare(d_cal,d_ref)
					d_real=[]
					for i in p1:
						d_real.append(d_ref[i])

					#比较后得到hkl(family)
					global hkl
					hkl=[]
					for i in p1:
						hkl.append(hkl_ref[i])

					######################################
					#from hkl family, get all hkl
					#to create HKL1，点三下中第一个点的等效族

					if A==B and A!=C and X==Y and float('%.4f'% X)==1.5708 and float('%.4f'% Z)==2.0944:
						#print('f*ck!!!!!!!!!!!!!!!')
						HKIL1=[hkl[0][0],hkl[0][1],(-1)*(hkl[0][0]+hkl[0][1]),hkl[0][2]]
						HKIL2=[hkl[1][0],hkl[1][1],(-1)*(hkl[1][0]+hkl[1][1]),hkl[1][2]]
						HKIL3=[hkl[2][0],hkl[2][1],(-1)*(hkl[2][0]+hkl[2][1]),hkl[2][2]]

						HKIL1=hkilfamily(HKIL1)
						HKIL2=hkilfamily(HKIL2)
						HKIL3=hkilfamily(HKIL3)

						HKL1=[]
						HKL2=[]
						HKL3=[]
							
						for i in range(len(HKIL1)):
							HKL1.append([HKIL1[i][0],HKIL1[i][1],HKIL1[i][3]])

						for i in range(len(HKIL2)):
							HKL2.append([HKIL2[i][0],HKIL2[i][1],HKIL2[i][3]])

						for i in range(len(HKIL3)):
							HKL3.append([HKIL3[i][0],HKIL3[i][1],HKIL3[i][3]])

						HKL1=delrepeat(HKL1)
						HKL2=delrepeat(HKL2)
						HKL3=delrepeat(HKL3)

					else:		 
						HKL1=hklfamily(hkl[0])
						if A==B:
							HKLAB=hklfamily([hkl[0][1],hkl[0][0],hkl[0][2]])
							HKL1=HKL1+HKLAB

						if B==C:
							HKLBC=hklfamily([hkl[0][0],hkl[0][2],hkl[0][1]])
							HKL1=HKL1+HKLBC

						if A==C:
							HKLAC=hklfamily([hkl[0][2],hkl[0][1],hkl[0][0]])
							HKL1=HKL1+HKLAC

						HKL1=delrepeat(HKL1)

						######to create HKL2，第二个点的等效族
						HKL2=hklfamily(hkl[1])
						if A==B:
							HKLAB=hklfamily([hkl[1][1],hkl[1][0],hkl[1][2]])
							HKL2=HKL2+HKLAB

						if B==C:
							HKLBC=hklfamily([hkl[1][0],hkl[1][2],hkl[1][1]])
							HKL2=HKL2+HKLBC

						if A==C:
							HKLAC=hklfamily([hkl[1][2],hkl[1][1],hkl[1][0]])
							HKL2=HKL2+HKLAC

						HKL2=delrepeat(HKL2)

						######to create HKL3，第三个点的等效族
						HKL3=hklfamily(hkl[2])
						if A==B:
							HKLAB=hklfamily([hkl[2][1],hkl[2][0],hkl[2][2]])
							HKL3=HKL3+HKLAB

						if B==C:
							HKLBC=hklfamily([hkl[2][0],hkl[2][2],hkl[2][1]])
							HKL3=HKL3+HKLBC

						if A==C:
							HKLAC=hklfamily([hkl[2][2],hkl[2][1],hkl[2][0]])
							HKL3=HKL3+HKLAC
						HKL3=delrepeat(HKL3)

					#############################################
					#有了点三下的三个点的等效晶面族，计算所有可能夹角
					#算点1，2
					angle_com1=[]
					t1=0

					#跑循环,所有跟另外的面的夹角
					while t1<len(HKL1):
						angle_com1.append([])
						m1=0
						while m1<len(HKL2):
						   p1=plane_angle(HKL1[t1],HKL2[m1],d_real[0],d_real[1])
						   angle_com1[t1].append(p1)
						   m1=m1+1
						t1=t1+1

					#这个list中：[[点1的可能],[点2的可能]]
					pre_direction1=find_close(angle_plot1,angle_com1)

					##################################################
					#确定当前晶带
					k4=0
					pre_zone1=[]
					while k4<len(pre_direction1[0]):
						z1=[HKL1[pre_direction1[0][k4]],HKL2[pre_direction1[1][k4]]]	#HKL1[最接近的项]，HKL2[最接近的项]
						pre_zone1.append(ortho(z1))									 #找垂直他们的晶带
						k4=k4+1

					pre_zone1=delpnzone(pre_zone1)			#去除重复

					for i in range(len(pre_zone1)):
						pre_zone1[i]=simplify(pre_zone1[i])	#化简

					#####################

					####################
					fin_zone1=[i for i in pre_zone1]			  	#存在
					fin_zone2=delrepeat(fin_zone1)					#去重，结果存于fin_zone2

					for i in range(len(fin_zone2)):
						if fin_zone2[i][0]>=0 and fin_zone2[i][1]>=0 and fin_zone2[i][2]>=0:
							fin_zone2=[fin_zone2[i]]
							break

					tem=[]
					tem.append(fin_zone2)

					for i in range(len(fin_zone2)):
						order_of_fin=fin_zone1.index(fin_zone2[i])
						OP_HKL=HKL1[pre_direction1[0][order_of_fin]]
						OQ_HKL=HKL2[pre_direction1[1][order_of_fin]]
						temangle=round(plane_angle(OP_HKL,OQ_HKL,d_real[0],d_real[1])*180/math.pi,2)
						tem.append(temangle)

					tem.append(round(angle_plot1*180/math.pi,2))

					print(tem)
					print(hkl)

					if abs(tem[-2]-tem[-1])>3.5:
						tem=[['ruled out'],round(angle_plot1*180/math.pi,2)]
						#print('killed')
					else:
						pass

					return tem
				####################################################################################
				#***************以上indexzone******************************************************#

				for yrw in range(5):   #以下5次循环跑5组点的标定结果
					global fin_zone
					fin_zone=[]
					HKL_list=[]

					three_points=[]
					three_points=[ten_points[2*yrw],ten_points[2*yrw+1],ten_points[2*yrw]]
						
					for dx in range(40):
						for dy in range(40):
							try:
								ttem=index_zone(three_points,(90+dx*0.5)/100,(90+dy*0.5)/100)
								
								for i in ttem[0]:
									fin_zone.append(i)							

							except Exception as e:
								fin_zone.append(str('err@@ ')+str(e))

							finally:
								fin_zone=delrepeat(fin_zone)

					print(str(yrw+1)+':  '+str(fin_zone))
					fin_yrw.append(fin_zone)

				print(fin_yrw)

				#####################################################################
				#用户选择五组标定结果###########################
				#********************************************#
				#回到状态2
				def dosame():
					global state
					state=2
					root.destroy()

				#回到状态1
				def doitagain():
					global state
					state=1
					root.destroy()
				
				#用户点击某个合格的按钮
				def choosepair(a):
					global state
					global yrww
					yrww=int(float(a))
					state=0
					root.destroy()
					
				root=tkinter.Tk()
				root.title('Result')
				root.geometry('800x550')

				#tkinter.Label(root,text='所选点的代表的晶向依次为 ').pack() 展示结果
				tkinter.Label(root,text='\nPlease choose the most probable current crystal zone:').pack()			 #展示结果

				tkinter.Button(root,text='First Pair:  '+str(fin_yrw[0]),command=lambda:choosepair(0)).pack()
				tkinter.Label(root,text='').pack()
				#tkinter.Label(root,text='indexd points'+str(showhkl[0]).replace('[','(').replace(']',')')).pack()
				#tkinter.Label(root,text='theory angles '+str(showangle[0][:-1])+' , measured angle='+str(showangle[0][-1])+'\n').pack()

				tkinter.Button(root,text='Second Pair:  '+str(fin_yrw[1]),command=lambda:choosepair(1)).pack()
				tkinter.Label(root,text='').pack()
				#tkinter.Label(root,text='indexd points'+str(showhkl[1]).replace('[','(').replace(']',')')).pack()
				#tkinter.Label(root,text='theory angles '+str(showangle[1][:-1])+' , measured angle='+str(showangle[1][-1])+'\n').pack()

				tkinter.Button(root,text='Third Pair:  '+str(fin_yrw[2]),command=lambda:choosepair(2)).pack()
				tkinter.Label(root,text='').pack()
				#tkinter.Label(root,text='indexd points'+str(showhkl[2]).replace('[','(').replace(']',')')).pack()
				#tkinter.Label(root,text='theory angles '+str(showangle[2][:-1])+' , measured angle='+str(showangle[2][-1])+'\n').pack()
				
				tkinter.Button(root,text='Forth Pair:  '+str(fin_yrw[3]),command=lambda:choosepair(3)).pack()
				tkinter.Label(root,text='').pack()
				#tkinter.Label(root,text='indexd points'+str(showhkl[3]).replace('[','(').replace(']',')')).pack()
				#tkinter.Label(root,text='theory angles '+str(showangle[3][:-1])+' , measured angle='+str(showangle[3][-1])+'\n').pack()

				tkinter.Button(root,text='Fifth Pair:  '+str(fin_yrw[4]),command=lambda:choosepair(4)).pack()
				tkinter.Label(root,text='').pack()
				#tkinter.Label(root,text='indexd points'+str(showhkl[4]).replace('[','(').replace(']',')')).pack()
				#tkinter.Label(root,text='theory angles '+str(showangle[4][:-1])+' , measured angle='+str(showangle[4][-1])).pack()
				
				tkinter.Label(root,text='').pack()
				tkinter.Button(root,text='Something is wrong. Select spots again.',command=dosame).pack()
				tkinter.Label(root,text='').pack()
				tkinter.Button(root,text='  Nope. Return to the beginning.  ',command=doitagain).pack()

				Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)
				tkinter.Label(root,text='').pack(side=tkinter.BOTTOM)
				Button(root,text='Any inconvience? Submit feedback',command=leave_comment).pack(side=tkinter.BOTTOM)

				Log.append('e')
				root.mainloop()	 #弹框

				#标定结束，弹窗关闭#############################
				#********************************************#
 
		#######################################################
		#标定完成，接下来是第二步：输入想要的晶带HKL，计算旋转角度

		three_points=[ten_points[2*yrww],ten_points[2*yrww+1],ten_points[2*yrww]]

		vital=[]

		for dx in range(40):
			vital.append([])
			for dy in range(40):
				try:
					ttem=index_zone(three_points,(90+dx*0.5)/100,(90+dy*0.5)/100)
					
					for i in ttem[0]:
						fin_zone.append(i)
						vital[dx].append(i)					

				except:
					fin_zone.append(str('err@'))
					vital[dx].append(str('err@'))

				finally:
					fin_zone=delrepeat(fin_zone)

		#print(vital)

		#保存标定结果
		with open("comment.txt","a") as f:
			f.write(str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))+'\n')
			f.write('Clicked on: '+str(three_points))
			f.write('Center: '+str(center))
			f.write('Indexed zone: ('+str(yrww+1)+'), out of '+str(fin_yrw)+' \n\n')

		#先确定OA，OB在图中的坐标

		def number_mul_vector(a,b):
			return [a*b[0],a*b[1],a*b[2]]
		def vector_add_vector(a,b):
			return [a[0]+b[0],a[1]+b[1],a[2]+b[2]]
		def vec_mul_vec(a,b):
			return a[0]*b[0]+a[1]*b[1]


		###############################修改旋转角正负

		OP_plot=[three_points[0][0]-center[0],three_points[0][1]-center[1]]	 ######
		OQ_plot=[three_points[1][0]-center[0],three_points[1][1]-center[1]]	 #将不太好的对称过来
						

#####################################################################################
		
		if len(fin_zone)==1:		   #此处是人手动输入识别出的多个晶带中正确的那个
			fin_zone4=fin_zone[0]
		else:
			def get_f4():
				global f4h
				global f4k
				global f4l
				f4h=etilt_f4_h.get()
				f4k=etilt_f4_k.get()
				f4l=etilt_f4_l.get()
				root.destroy()
						
			root=tkinter.Tk()
			root.title('Choose the right one')
			root.geometry('700x360')
			tkinter.Label(root,text='').pack()
			tkinter.Label(root,text='Possible zones are \n'+str(fin_zone)).pack()
			tkinter.Label(root,text='').pack()

			def choosef3(a):
				global YRW
				YRW=int(float(a))
				root.destroy()

			for i in range(len(fin_zone)):#循环显示button
				Button(root,text='It is '+str(fin_zone[i]),command=lambda i=i:choosef3(i)).pack()
			

			Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)

			Log.append('f')
			root.mainloop()

			fin_zone4=fin_zone[YRW]

			dx=0
			while dx<40:
				dy=0
				while dy<40:
					if vital[dx][dy]==fin_zone4:
						ttem=index_zone(three_points,(90+dx*0.5)/100,(90+dy*0.5)/100)
						print('!  '+str(ttem[0]))

						dx=100
						dy=100
					else:
						dy+=1
				dx+=1


		
		while state==0:
			#弹窗输入目前双倾杆位置，放于x1y1##########
			#**************************************#
			root=tkinter.Tk()
			root.title('Input current holder position')
			root.geometry('500x200')
			tkinter.Label(root,text='\nPlease give the current tilt angles (x&y) of the holder.(degree, 0-180)\n\n\n').pack()

			etiltx=Entry(root,width=5)
			etiltx.place(x=200,y=43)
			etilty=Entry(root,width=5)
			etilty.place(x=270,y=43)

			def holder_position():
				global x1
				global y1
				x1=etiltx.get()
				y1=etilty.get()
				x1=float(x1)/180*math.pi
				y1=float(y1)/180*math.pi
				root.destroy()

			button=tkinter.Button(root,text='Done!',command=holder_position)
			button.pack()

			Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)
			root.mainloop()

			#**************************************#
			#弹窗结束################################

			###############################
			#UI结束，所有信息录入。开始计算带轴到带轴，用推导的经-纬模型下的公式
			###############################

			############################################
			def number_mul_vector(a,b):
				return [a*b[0],a*b[1],a*b[2]]
			def vector_add_vector(a,b):
				return [a[0]+b[0],a[1]+b[1],a[2]+b[2]]
			def vec_mul_vec(a,b):
				return a[0]*b[0]+a[1]*b[1]
			############################################
			############################################
			
			def zone_to_zone(tar_zone):
				#calculate theta
				#确定OP，OQ所代表的HKL
				global OP_HKL
				global OQ_HKL
				global fin_zone1
				order_of_fin=fin_zone1.index(fin_zone4)
				OP_HKL=HKL1[pre_direction1[0][order_of_fin]]
				OQ_HKL=HKL2[pre_direction1[1][order_of_fin]]

				#计算OA，OB,OF所代表的HKL
				angle_POA=measureangle(center,three_points[0],point_A)
				angle_QOA=measureangle(center,three_points[1],point_A)
				angle_POB=measureangle(center,three_points[0],point_B)
				angle_QOB=measureangle(center,three_points[1],point_B)
				angle_POQ=measureangle(center,three_points[0],three_points[1])
				OP_length=distance(three_points[0],center)
				OQ_length=distance(three_points[1],center)
				OF_HKL=ortho([tar_zone,fin_zone4])

				#将A、B投影在OP和OQ基矢量上（分解）
				OA_HKL_1=number_mul_vector(abs(50*(2)**0.5/OP_length*math.sin(angle_QOA)/math.sin(angle_POQ)),OP_HKL)
				OA_HKL_2=number_mul_vector(abs(50*(2)**0.5/OQ_length*math.sin(angle_POA)/math.sin(angle_POQ)),OQ_HKL)
				OB_HKL_1=number_mul_vector(abs(50*(2)**0.5/OP_length*math.sin(angle_QOB)/math.sin(angle_POQ)),OP_HKL)
				OB_HKL_2=number_mul_vector(abs(50*(2)**0.5/OQ_length*math.sin(angle_POB)/math.sin(angle_POQ)),OQ_HKL)
				
				global OA_HKL
				global OB_HKL
				OA_HKL=[]
				OB_HKL=[]

				if vec_mul_vec(OA_plot,OP_plot)>0:
					pass
				else:
					OA_HKL_1=number_mul_vector(-1,OA_HKL_1)
				if vec_mul_vec(OA_plot,OQ_plot)>0:
					pass
				else:
					OA_HKL_2=number_mul_vector(-1,OA_HKL_2)
				OA_HKL=vector_add_vector(OA_HKL_1,OA_HKL_2)

				if vec_mul_vec(OB_plot,OP_plot)>0:
					pass
				else:
					OB_HKL_1=number_mul_vector(-1,OB_HKL_1)
				if vec_mul_vec(OB_plot,OQ_plot)>0:
					pass
				else:
					OB_HKL_2=number_mul_vector(-1,OB_HKL_2)
				OB_HKL=vector_add_vector(OB_HKL_1,OB_HKL_2)


				#计算A，B，F点所代表的d
				d_point_A=((OA_HKL[0]/A/math.sin(Y))**2+(OA_HKL[2]/C/math.sin(Y))**2+2*OA_HKL[0]*OA_HKL[1]*math.cos(Y)/A/C/math.sin(Y)**2+OA_HKL[1]**2/B**2)**(-0.5)
				d_point_F=((OF_HKL[0]/A/math.sin(Y))**2+(OF_HKL[2]/C/math.sin(Y))**2+2*OF_HKL[0]*OF_HKL[1]*math.cos(Y)/A/C/math.sin(Y)**2+OF_HKL[1]**2/B**2)**(-0.5)
				d_point_B=((OB_HKL[0]/A/math.sin(Y))**2+(OB_HKL[2]/C/math.sin(Y))**2+2*OB_HKL[0]*OB_HKL[1]*math.cos(Y)/A/C/math.sin(Y)**2+OB_HKL[1]**2/B**2)**(-0.5)			

				#利用 文献中的由旋转矩阵推得的结果 计算杆的旋转角度 
				global MON_angle
				global FOA_angle
				MON_angle=zone_angle(fin_zone4,tar_zone)  #这里应该是晶向夹角
				FOA_angle=plane_angle(OF_HKL,OA_HKL,d_point_F,d_point_A)

				normal_x=math.sin(MON_angle)*math.sin(FOA_angle)
				normal_y=math.sin(MON_angle)*math.cos(FOA_angle)
				normal_z=math.cos(MON_angle)
				normal_w=1
				
				OA_rotate=-(arcsin(math.sin(MON_angle)*math.cos(math.pi-FOA_angle)*math.cos(x1)+math.sin(x1)*math.cos(MON_angle))-x1)
				OB_rotate=-(arctan(tan(MON_angle)*math.sin(math.pi-FOA_angle)/(math.cos(x1)-math.sin(x1)*tan(MON_angle)*math.cos(math.pi-FOA_angle))))

				global OA_rotate_new
				global OB_rotate_new
				global OA_rotate_new2
				OA_rotate_new=arcsin((normal_z*math.sin(x1)+normal_y*math.cos(x1))/normal_w)-x1
				OB_rotate_new=-arctan(normal_x*math.cos(x1)/(normal_z-(normal_z*math.sin(x1))**2-normal_z*normal_y*math.sin(x1)*math.cos(x1)))

				OA_rotate_new2=(normal_z*math.sin(x1)+normal_y*math.cos(x1))/(normal_z*math.cos(x1)-normal_y*math.sin(x1))*math.cos(-OB_rotate_new)-x1

				global x2_sol
				global y2_sol
				x2_sol=OA_rotate
				y2_sol=OB_rotate

				def theta(m,n):#from the rotate angle to the tilt angle,m=x2,n=y2 
					return math.acos(math.cos(x1)*math.cos(m)*math.cos(y1)*math.cos(n)+math.cos(x1)*math.cos(m)*math.sin(y1)*math.sin(n)+math.sin(x1)*math.sin(m))
	 
				global pre_MON
				pre_MON=theta(x2_sol+x1,y2_sol+y1)#目前的旋转角

			#################################################################
			
			#自定义zone输入弹窗######################
			#弹窗输入自定义位置，放于tar_zone#########
			#**************************************#
			showbasicX=[]
			showbasicX2=[]
			showbasicY=[]
			basicset=[[0,0,1],[0,1,0],[1,0,0],[0,1,1],[1,0,1],[1,1,0],[1,1,1]]
			for i in basicset:
				try:
					zone_to_zone(i)
					showbasicX.append(round((OA_rotate_new+x1)/math.pi*180,2))
					showbasicX2.append(round((OA_rotate_new2+x1)/math.pi*180,2))
					showbasicY.append(round((OB_rotate_new+y1)/math.pi*180,2))
				except:
					showbasicX.append('err')
					showbasicX2.append('err')
					showbasicY.append('err')

			root=tkinter.Tk()
			root.title('Next')
			root.geometry('700x400')
			tkinter.Label(root,text='\nResults of defult zones: '+str(basicset)).pack()
			tkinter.Label(root,text='\n Final position of x: '+str(showbasicX)).pack()
			#tkinter.Label(root,text='\n(arctan)Final position of x: '+str(showbasicX2)).pack()
			tkinter.Label(root,text='\n	Final position of y: '+str(showbasicY)).pack()

			tkinter.Label(root,text='\nOr give the desired zone').pack()
			tkinter.Label(root,text='\n\n[	,	,	]\n').pack()

			Hconst=StringVar()
			Hconst.set('0')		   #缺省值
			etarH=Entry(root,textvariable=Hconst,width=5)
			etarH.place(x=280,y=200)

			Kconst=StringVar()
			Kconst.set('0')		   #缺省值			
			etarK=Entry(root,textvariable=Kconst,width=5)
			etarK.place(x=330,y=200)

			Lconst=StringVar()
			Lconst.set('1')		   #缺省值
			etarL=Entry(root,textvariable=Lconst,width=5)
			etarL.place(x=380,y=200)

			def get_target_zone():
				global tar_zone	

				tar_zoneH=etarH.get()
				tar_zoneK=etarK.get()
				tar_zoneL=etarL.get()
				tar_zone=[]
				tar_zone.append(int(float(tar_zoneH)))
				tar_zone.append(int(float(tar_zoneK)))
				tar_zone.append(int(float(tar_zoneL)))
				
				zone_to_zone(tar_zone)

				window=tkinter.Tk()
				window.title('Final Results')
				window.geometry('700x300')

				def windowclose():
					window.destroy()

				Label(window,text='\nThe current zone is '+str(fin_zone4)+' , the target zone is '+str(tar_zone)).pack()
				Label(window,text='The current tilt angles of double-tilt holder').pack()
				Label(window,text='x1=  '+str(round(x1/math.pi*180,2))+' ,   y1=  '+str(round(y1/math.pi*180,2))).pack()

				#Label(window,text='\nBased on the given information, the rotating angles of x& y are').pack()
				#Label(window,text='△x=  '+str(round(x2_sol/math.pi*180,2))+' ,   △y=  '+str(round(y2_sol/math.pi*180,2))).pack()
				Label(window,text='Total rotating angle = '+str(round(pre_MON/math.pi*180,2))+'. The angle between the two zones = '+str(round(MON_angle/math.pi*180,2))).pack()
			
				#Label(window,text='\nFinal position of the holder x2=  '+str(round((x2_sol+x1)/math.pi*180,2))+' , y2=  '+str(round((y2_sol+y1)/math.pi*180,2))).pack()
				Label(window,text='\n Final position of the holder x2 =  '+str(round((OA_rotate_new+x1)/math.pi*180,2))+' , y2=  '+str(round((OB_rotate_new+y1)/math.pi*180,2))).pack()
				#Label(window,text='\n(arctan)Final position of the holder x2 =  '+str(round((OA_rotate_new2+x1)/math.pi*180,2))+' , y2=  '+str(round((OB_rotate_new+y1)/math.pi*180,2))).pack()
				Label(window,text='').pack()

				Button(window,text='Well done, close window',command=windowclose).pack(side=tkinter.BOTTOM)
				
				window.mainloop()
				
			button=tkinter.Button(root,text='Done!',command=get_target_zone)
			button.pack()
			Label(root,text='').pack()

			def closemain():
				root.destroy()

			button=tkinter.Button(root,text='NEXT',command=closemain)
			button.pack()

			Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)

			Log.append('g')
			root.mainloop()
			#**************************************#
			###################################################################
			
			#弹窗显示答案####################
			#******************************#
			root=tkinter.Tk()
			root.title('Final')
			root.geometry('500x250')

			#state=0，后半循环继续，回到while state=0
			def backtof1():
				global state
				state=0
				root.destroy()

			#state≠0，后半循环结束，回到开头while true	
			def backtof2():
				global state
				state=1
				root.destroy()

			####################################
			Label(root,text='').pack()
			Button(root,text='Nope，wrong input of the holder position',command=backtof1).pack()
			
			Label(root,text='').pack(side=tkinter.BOTTOM)			
			Button(root,text='Any inconvience? Submit feedback',command=leave_comment).pack(side=tkinter.BOTTOM)
			Label(root,text='').pack(side=tkinter.BOTTOM)
			Button(root,text='Return to beginning',command=backtof2).pack(side=tkinter.BOTTOM)
			Label(root,text='').pack(side=tkinter.BOTTOM)
			Button(root,text='Well done, quit',command=goodbye).pack(side=tkinter.BOTTOM)

			'''
			with open("comment.txt","a") as f:  #日志记录结果
				f.write(str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))+'\n')
				f.write('OA_HKL: '+str(OA_HKL)+' OB_HKL: '+str(OB_HKL)+' OP_HKL: '+str(OP_HKL)+' OQ_HKL: '+str(OQ_HKL)+'\n')
				f.write('OP_plot: '+str(OP_plot)+' OQ_plot: '+str(OQ_plot)+'\n')
				f.write('Current zone: '+str(fin_zone4)+' Target zone: '+str(tar_zone)+'\n')
				f.write('Tilt:	x1= '+str(round(x1/math.pi*180,2))+' ,  y1= '+str(round(y1/math.pi*180,2))+'\n')
				f.write('Result: △x= '+str(round(x2_sol/math.pi*180,2))+' , △y=  '+str(round(y2_sol/math.pi*180,2))+'\n\n')
			
			'''

			Log.append('h')
			root.mainloop()
			#******************************#
			#弹窗关闭########################
			
			###########################################
			#final

###########################
###异常情况弹窗
	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		from tkinter import *
		def goodbye():
			with open("comment.txt","a") as f:
				f.write('LOG: '+str(Log)+'\n')
				f.write('*********************\n')
			sys.exit(0)

		with open("comment.txt","a") as f:	  #保存错误记录
			f.write(str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))+'\n')
			f.write('Error: line '+str(exc_tb.tb_lineno)+', '+str(exc_type)+', '+str(e)+'\n\n')
		
		root=tkinter.Tk()
		root.title('Ooops!')
		root.geometry('500x200')
		Label(root,text='\n Error！').pack()
		Label(root,text='\n'+str(e)+'\n').pack()
		Button(root,text='Fine',command=root.destroy).pack()
		Button(root,text='emergency exit',command=goodbye).pack(side=tkinter.BOTTOM)

		Log.append('Y')
		root.mainloop()
		
