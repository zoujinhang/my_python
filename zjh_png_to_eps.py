'''
新系统需要设置convert权限
sudo gedit /etc/ImageMagick-6/policy.xml
添加或修改如下
<policy domain="coder" rights="read|write" pattern="EPS" />
<policy domain="coder" rights="read|write" pattern="PDF" />
<policy domain="coder" rights="read|write" pattern="XPS" />
<policy domain="coder" rights="read|write" pattern="PNG" />
'''
import os

def png_to_eps(fit,out):
	cmd = 'convert -density 3000 '+fit+' '+out
	os.system(cmd)
#fit = '/home/laojin/result/bn200219317_zbb/B_txx_15-350.png'
#out = '/home/laojin/result/bn200219317_zbb/B_txx_15-350.eps'
#png_to_eps(fit,out)

