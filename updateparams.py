def setparamsfile_CLASS (params_list, pfilename, pnewfilename):
	with open (pfilename) as file: lines = file.read ().splitlines ()
	with open (pnewfilename, 'w') as file:
		for line in lines:
			if line.startswith ('h '):
				file.write (line.replace (line, 'h = {0}' .format (params_list[0])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[0]))
			
			elif line.startswith ('Omega_b'):
				file.write (line.replace (line, 'Omega_b = {0}' .format (params_list[1])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[1]))
			
			elif line.startswith ('Omega_cdm'):
				file.write (line.replace (line, 'Omega_cdm = {0}' .format (params_list[2])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[2]))
			
			elif line.startswith ('Omega_Lambda'):
				file.write (line.replace (line, 'Omega_Lambda = {0}' .format (params_list[3])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[3]))
			
			elif line.startswith ('tau_reio'):
				file.write (line.replace (line, 'tau_reio = {0}' .format (params_list[4])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[4]))
			
			elif line.startswith ('A_s'):
				file.write (line.replace (line, 'A_s = {0}' .format (params_list[5])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[5]))
			
			elif line.startswith ('n_s'):
				file.write (line.replace (line, 'n_s = {0}' .format (params_list[6])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[6]))
			
			elif line.startswith ('N_ur'):
				file.write (line.replace (line, 'N_ur = {0}' .format (params_list[7])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[7]))
			
			elif line.startswith ('N_cdm'):
				file.write (line.replace (line, 'N_cdm = {0}' .format (params_list[8])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[8]))
			
			elif line.startswith ('m_cdm'):
				file.write (line.replace (line, 'm_cdm = {0}' .format (params_list[9])))
				file.write ('\n')
				print (line, '-> {0}'.format (params_list[9]))
			
			elif line.startswith ('T_ncdm'):
				if not params_list[6] == 0:
					file.write (line.replace (line, 'T_ncdm = {0}' .format (params_list[6])))
					file.write ('\n')
					print (line, '-> {0}'.format (params_list[6]))
			
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
				new = params_list[1]*params_list[0]**2
				file.write (lines[i].replace (lines[i], '{0}' .format (new)))
				file.write ('\n')
				print (old, '-> {0}'.format (new))
			elif i == 2:
				old = lines[i]
				new = (params_list[1]+params_list[2])*params_list[0]**2
				file.write (lines[i].replace (lines[i], '{0}' .format (new)))
				file.write ('\n')
				print (old, '-> {0}'.format (new))
			elif i == 4:
				old = lines[i]
				new = params_list[3]*params_list[0]**2
				file.write (lines[i].replace (lines[i], '{0}' .format (new)))
				file.write ('\n')
				print (old, '-> {0}'.format (new))
			elif i == 7:
				if not params_list[7] == 3.046:
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
