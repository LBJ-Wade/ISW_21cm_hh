def setparamsfile_CLASS (params_list, pfilename, pnewfilename):
	with open (pfilename) as file: lines = file.read ().splitlines ()
	with open (pnewfilename, 'w') as file:
		for line in lines:
			if line.startswith ('omega_cdm'):
				file.write (line.replace (line, 'omega_cdm = {0}' .format (params_list[0])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[0]))
			
			elif line.startswith ('omega_b'):
				file.write (line.replace (line, 'omega_b = {0}' .format (params_list[1])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[1]))
			
			elif line.startswith ('100*theta_s'):
				file.write (line.replace (line, '100*theta_s = {0}' .format (params_list[2])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[2]))
			
			elif line.startswith ('tau_reio'):
				file.write (line.replace (line, 'tau_reio = {0}' .format (params_list[3])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[3]))
			
			elif line.startswith ('A_s'):
				file.write (line.replace (line, 'A_s = {0}' .format (params_list[4])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[4]))
			
			elif line.startswith ('n_s'):
				file.write (line.replace (line, 'n_s = {0}' .format (params_list[5])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[5]))
			
			elif line.startswith ('m_ncdm'):
				file.write (line.replace (line, 'm_ncdm = {0}' .format (params_list[6])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[6]))
			
			elif line.startswith ('T_ncdm'):
				file.write (line.replace (line, 'T_ncdm = {0}' .format (params_list[7])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[7]))
			
			else:
				print (line)
				file.write(line)
				file.write('\n')

def setparamsfile_HYREC (params_list, pfilename, pnewfilename):
	with open (pfilename) as file: lines = file.read ().splitlines ()
	with open (pnewfilename, 'w') as file:
		for i in range(len(lines)):
			if i == 1:
				old = lines[i]
				new = params_list[1]
				file.write (lines[i].replace (lines[i], '{0}' .format (new)))
				file.write ('\n')
				print (old, '-> {0}'.format (new))
			elif i == 2:
				old = lines[i]
				new = (params_list[0]+params_list[1])
				file.write (lines[i].replace (lines[i], '{0}' .format (new)))
				file.write ('\n')
				print (old, '-> {0}'.format (new))
			elif i == 4:
				old = lines[i]
				new = params_list[8]
				file.write (lines[i].replace (lines[i], '{0}' .format (new)))
				file.write ('\n')
				print (old, '-> {0}'.format (new))
			elif i == 7:
				if not params_list[9] == 3.046:
					old = lines[i]
					new = params_list[7]
					file.write (lines[i].replace (lines[i], '{0}' .format (new)))
					file.write ('\n')
					print (old, '-> {0}'.format (new))
				else:
					print (lines[i])	
			else:
				print (lines[i])
				file.write(lines[i])
				file.write('\n')
