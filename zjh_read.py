from astropy.io import fits

class Read_tte_file(object):
	def __init__(self,dailydata = None,triggerdata = None):
		if dailydata is not None :
			self.file = fits.open(dailydata)
			self.time = self.file[2].data.field(0)
			self.ch = self.file[2].data.field(1)
			self.ch_n = self.file[1].data.field(0)
			self.e1 = self.file[1].data.field(1)
			self.e2 = self.file[1].data.field(2)
		if triggerdata is not None:
			self.file = fits.open(triggerdata)
			self.trigtime = self.file[0].header['TRIGTIME']
			self.time = self.file[2].data.field(0)
			self.ch = self.file[2].data.field(1)
			self.ch_n = self.file[1].data.field(0)
			self.e1 = self.file[1].data.field(1)
			self.e2 = self.file[1].data.field(2)
	def get_file(self):

		return self.file

	def get_time_ch(self):

		return self.time,self.ch

	def get_ch_n_e1_e2(self):

		return self.ch_n,self.e1,self.e2


















