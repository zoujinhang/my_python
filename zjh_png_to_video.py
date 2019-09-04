import cv2

def png_to_video(linklist,save,fps,size,topdir = ''):
    video = cv2.VideoWriter(topdir+save,cv2.VideoWriter_fourcc(*'XVID'),fps,size)
    for item in linklist:
        print('------------------------')
        print(item)
        img = cv2.imread(topdir+item)
        video.write(img)
        print('------------------------')
    video.release()












