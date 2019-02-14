

def anyOverlap(s1, e1, s2, e2):

	"""
	Determine whether or not the interval s1-e1 has any overlap with
	the interval s2-e2
	"""
	
	ovrlp = False
	
	s1 = int(s1)
	e1 = int(e1)
	s2 = int(s2)
	e2 = int(e2)

	# Return 'no overlap' unless:

	# s1 OR e1 is between s2-e2
	if (s2 <= s1 <= e2) | (s2 <= e1 <= e2):
		ovrlp = True

	# if s2-e2 is inside s1-e1
	elif (s1 <= s2) & (e1 >= e2):
		ovrlp = True
	
	return(ovrlp)




def createTranscriptMapper(exons, sCDS, eCDS, strand):

	"""
	Create a MAPPER for the codons of a transcript. It is a dictionary with the
	positions and codons of every exon in the transcript.
	"""

	mapper = {}

	# set a variable that indicates the status of the CDS
	cds = 'not_started'

	# loop over exons
	for i in range(len(exons)):

		# get exon number and create a subdictionary for the exon
		EXON = i + 1
		mapper[EXON] = {}

		# set START and END position of the exon
		sE = int(exons[i][0])
		eE = int(exons[i][1])
		mapper[EXON]['position'] = (sE, eE)

		# set as zero the other info for this exon
		mapper[EXON]['codon_position'] = (0, 0)
		mapper[EXON]['codon_number'] = (0, 0)
		mapper[EXON]['extra_bases'] = (0, 0)

		# if the the previous exons did not contain CDS
		if cds == 'not_started':

			# if the current exon does not contain CDS as well
			if sCDS > eE:
				continue

			# if there is a CDS
			else:

				# change CDS status
				cds = 'started'

				# count codons between start of CDS and end of exon
				l = eE - sCDS
				n_cod = l//3

				# in case of extra bases at 3', add one codon
				extraEnd = l%3
				if extraEnd > 0: 
					n_cod = n_cod + 1

				# correct exon info
				mapper[EXON]['codon_position'] = (sCDS, eE)
				mapper[EXON]['codon_number'] = (1, n_cod)
				mapper[EXON]['extra_bases'] = (0, extraEnd)


		# if the CDS started in the previous exons and is not yet finished
		elif cds == 'started':

			# compute possible extra bases at 5'
			if extraEnd > 0: extraStart = 3 - extraEnd
			else: extraStart = 0

			### derive (1)codon positions and (2)size to end of CDS ###
			# (they change dipending on whether the CDS finishes in 
			# current exon)
			if eCDS > eE:
				l = eE - sE - extraStart
				mapper[EXON]['codon_position'] = (sE, eE)
			else:
				l = eCDS - sE - extraStart
				mapper[EXON]['codon_position'] = (sE, eCDS)
				
				# change CDS status if it finishes in this exon
				cds = 'ended'

			# start codon is the same as previous exon's end codon
			sCOD = mapper[EXON - 1]['codon_number'][1]

			# count number of codons until end of CDS
			n_cod = l//3
			eCOD = sCOD + n_cod

			# add 1 COD to start COD, in case there are NO extra bases at 5'
			if extraStart == 0: sCOD = sCOD + 1

			# compute extra bases at 3'
			extraEnd = l%3
			# add 1 COD if there are extra bases
			if extraEnd > 0: eCOD = eCOD + 1

			# correct exon info
			mapper[EXON]['codon_number'] = (sCOD, eCOD)
			mapper[EXON]['extra_bases'] = (extraStart, extraEnd)


		# if the CDS ended in previous exons, continue with no changes
		elif cds == 'ended':
			continue


	#### ADJUST THE MAPPER IN CASE OF NEGATIVE STRAND ####
	# (re-compute codon numbers in the opposite order (last codon is 1)
	if strand == '-':
		eCOD = 1
		mapper2 = {}

		# get reverse list of exon numbers
		keys = [i for i in mapper][::-1]

		# loop over exons in reverse
		for k in keys:

			# get reverse number for exon 
			k2 = len(mapper) - k + 1

			# copy info of current exon
			mapper2[k2] = mapper[k]

			# if the exon contains a CDS
			cod = mapper[k]['codon_number']
			if cod != (0, 0):
				# count the CODs and add them to the last COD of previous exon
				sCOD = cod[1] - cod[0] + eCOD
				# correct codon numbers
				mapper2[k2]['codon_number'] = (sCOD, eCOD)
				# set start COD as last COD to add to next exon				
				eCOD = sCOD
				# if this is not the last codon
				if k > 1:
					# add 1 COD in case end COD is shared with the next exon
					pre_cod = mapper[k - 1]['codon_number']
					if pre_cod[1] != cod[0]:
						eCOD = eCOD + 1

		mapper = mapper2

	return(mapper)



def computeCodons(pos, mapper, EXON, strand):

	sCpos, eCpos = mapper[EXON]['codon_position']
	sC, eC = mapper[EXON]['codon_number']
	sExtra, eExtra = mapper[EXON]['extra_bases']

	if sExtra > 0:
		sCpos = sCpos + sExtra
		if strand == '+':
			sC = sC + 1
		else:
			sC = sC - 1

	if eExtra > 0: 
		eCpos = eCpos - eExtra
		if strand == '+':
			eC = eC - 1
		else:
			eC = eC + 1


	UPl = pos - sCpos
	DOWNl = eCpos - pos

	UPcod = UPl//3
	UPcds = [sC, 0]

	if strand == '+':
		UPcds[1] = sC + UPcod
		if UPcod > 0:
			UPcds[1] = UPcds[1] - 1
	else:
		UPcds[1] = sC - UPcod
		if UPcod > 0:
			UPcds[1] = UPcds[1] + 1


	DOWNcod = DOWNl//3
	DOWNcds = [0, eC]

	if strand == '+':
		DOWNcds[0] = eC - DOWNcod
		if DOWNcod > 0:
			DOWNcds[0] = DOWNcds[0] + 1
	else:
		DOWNcds[0] = eC + DOWNcod
		if DOWNcod > 0:
			DOWNcds[0] = DOWNcds[0] - 1

	return((UPcds, DOWNcds))





def annotateSymbolTranscriptExonsCodons(bed, ucsc):
	
	import pandas as pd

	"""
	Annotate a BED file (standard, 3 columns) with HUGO symbol, every RefSeq
	transcript ID, every exon number and codons covered. Use a UCSC RefSeq
	reference table.
	"""

	cols = ['chrom', 'start', 'end', 'hugo', 'strand', 'refseq_transcript', 'exon', 'intron5', 'utr5', 'extra_bases5', 'first_complete_codon', 'last_complete_codon', 'extra_bases3', 'utr3', 'intron3', 'exon_start', 'cds_start', 'cds_end', 'exon_end']

	DF = pd.DataFrame(columns=cols)

	# loop over BED intervals (rows)
	for i, row in bed.iterrows():
		
		if i%10 == 0: print('%s/%s' %(i, len(bed)))

		# get coordinates of QUERY interval
		cQ = row.chrom
		sQ = row.start
		eQ = row.end

		# subset REF table for current chromosome
		ucscCHR = ucsc.loc[ucsc.chrom == cQ,].copy()

		# loop over genes/transcripts (rows) in chromosome
		for i2, row2 in ucscCHR.iterrows():

			# get transcript START and END
			sT = row2.txStart
			eT = row2.txEnd

			# if QUERY overlaps with the transcript
			if anyOverlap(sQ, eQ, sT, eT):
				
				df = pd.DataFrame(columns=cols)

				# get hugo and transcript ID and append to growing list
				hugo = row2.name2
				refseq = row2['name']
				strand = row2.strand

				# get every exon START and END position as list of tuples
				starts = row2.exonStarts.split(',')[:-1]
				ends = row2.exonEnds.split(',')[:-1]
				exons = list(zip(starts, ends))
				
				sCDS = int(row2.cdsStart)
				eCDS = int(row2.cdsEnd)


				# create MAPPER for transcript (dict to summarize codon info)
				mapper = createTranscriptMapper(exons, sCDS, eCDS, strand)

				# track when CDS starts and ends
				cds_status = 'pre_cds'

				# loop over exons:
				for EXON in mapper:

					sE, eE = mapper[EXON]['position']

					if anyOverlap(sQ, eQ, sE, eE):
						
						sCpos, eCpos = mapper[EXON]['codon_position']
						sC, eC = mapper[EXON]['codon_number']
						sExtra, eExtra = mapper[EXON]['extra_bases']

						if sExtra > 0: sC = sC + 1
						if eExtra > 0: eC = eC - 1

						if (sCpos != 0) & (sE < sCpos):
							cds_status = 'cds'

						if (eCpos != 0) & (eCpos < eE):
							cds_status = 'post_cds'


						intron5 = '0'
						utr5 = '0'
						extra5 = '0'
						cds = ['0', '0']
						extra3 = '0'
						utr3 = '0'
						intron3 = '0'


						if (sQ < sE) & (eQ > eE): 

							intron5 = str(sE - sQ)
							intron3 = str(eQ - eE)

							if (sCpos, eCpos) == (0, 0):
								if cds_status == 'pre_cds': utr5 = str(eE - sE)
								elif cds_status == 'post_cds': utr3 = str(eE - sE)

							else:
								if sE < sCpos:
									utr5 = str(sCpos - sE)
									cds = computeCodons(sCpos, mapper, EXON, strand)[1]
									extra3 = str(eExtra)

								elif eCpos < eE:
									extra5 = str(sExtra)
									cds = computeCodons(eCpos, mapper, EXON, strand)[0]
									utr3 = str(eE - eCpos)
									
								else:
									extra5 = str(sExtra)
									cds = (sC, eC)
									extra3 = str(eExtra)





						if (sQ < sE) & (eQ <= eE): 
							intron5 = str(sE - sQ)

							if (sCpos, eCpos) == (0, 0):
								if cds_status == 'pre_cds': utr5 = str(eQ - sE)
								elif cds_status == 'post_cds': utr3 = str(eQ - sE)
							else:
								if sE < sCpos:

									# compute extra3 bases in case eE is less
									# than 2 bases bigger than eQ
									if (eE - eQ < 2) & (eExtra > 0):
										extra3 = eExtra - eE + eQ

									if eQ <= sCpos:
										utr5 = str(eQ - sE)
									else:
										utr5 = str(sCpos - sE)
										cds = computeCodons(eQ, mapper, EXON, strand)[0]

								elif eCpos < eE:
									extra5 = str(sExtra)
									if eQ <= eCpos:
										cds = computeCodons(eQ, mapper, EXON, strand)[0]
									else:
										cds = computeCodons(eCpos, mapper, EXON, strand)[0]
										utr3 = str(eQ - eCpos)

								else:
									extra5 = str(sExtra)
									cds = computeCodons(eQ, mapper, EXON, strand)[0]
									# compute extra3 bases in case eE is less
									# than 2 bases bigger than eQ
									if (eE - eQ < 2) & (eExtra > 0):
										extra3 = eExtra - eE + eQ




						if (sQ >= sE) & (eQ > eE): 
							intron3 = str(eQ - eE)

#							if refseq == 'NM_000546.5':
#								if EXON == 3:
#									print('yes')


							if (sCpos, eCpos) == (0, 0):
								if cds_status == 'pre_cds': utr5 = str(eE - sQ)
								elif cds_status == 'post_cds': utr3 = str(eE - sQ)
							else:
								if sE < sCpos:
									extra3 = str(eExtra)
									if sQ < sCpos:
										utr5 = str(sCpos - sQ)
										cds = computeCodons(sCpos, mapper, EXON, strand)[1]
									else:
										cds = computeCodons(sQ, mapper, EXON, strand)[1]

								elif eCpos < eE:

									# compute extra5 bases in case sQ is less
									# than 2 bases bigger than sE
									if (sQ - sE < 2) & (sExtra > 0):
										extra5 = sExtra - sQ + sE
 
									if sQ < eCpos:
										cds = computeCodons(sQ, mapper, EXON, strand)[1]
										utr3 = str(eE - eCpos)
									else:
										utr3 = str(eE - sQ)

								else:
									# compute extra5 bases in case sQ is less
									# than 2 bases bigger than sE
									if (sQ - sE < 2) & (sExtra > 0):
										extra5 = sExtra - sQ + sE

									cds = computeCodons(sQ, mapper, EXON, strand)[1]
									extra3 = str(eExtra)
									



						if (sQ >= sE) & (eQ <= eE): 

							if (sCpos, eCpos) == (0, 0):
								if cds_status == 'pre_cds': utr5 = str(eQ - sQ)
								elif cds_status == 'post_cds': utr3 = str(eQ - sQ)

							else:
								if sE < sCpos:
									if eQ < sCpos:
										utr5 = str(eQ - sQ)
									elif sQ < sCpos < eQ:
										utr5 = str(sCpos - sQ)
										cds = computeCodons(eQ, mapper, EXON, strand)[0]
										# compute extra3 bases in case eE is
										# less than 2 bases bigger than eQ
										if (eE - eQ < 2) & (eExtra > 0):
											extra3 = eExtra - eE + eQ

									elif sCpos < sQ:
										cds[0] = computeCodons(sQ, mapper, EXON, strand)[1][0]
										cds[1] = computeCodons(eQ, mapper, EXON, strand)[0][1]
										# compute extra3 bases in case eE is
										# less than 2 bases bigger than eQ
										if (eE - eQ < 2) & (eExtra > 0):
											extra3 = eExtra - eE + eQ

								elif eCpos < eE:
									if eQ < eCpos:
										cds[0] = computeCodons(sQ, mapper, EXON, strand)[1][0]
										cds[1] = computeCodons(eQ, mapper, EXON, strand)[0][1]
										# compute extra5 bases in case sQ is
										# less than 2 bases bigger than sE
										if (sQ - sE < 2) & (sExtra > 0):
											extra5 = sExtra - sQ + sE

									elif sQ < eCpos < eQ:
										# compute extra5 bases in case sQ is
										# less than 2 bases bigger than sE
										if (sQ - sE < 2) & (sExtra > 0):
											extra5 = sExtra - sQ + sE
										cds = computeCodons(sQ, mapper, EXON, strand)[1]
										utr3 = str(eQ - eCpos)
									elif eCpos < sQ:
										utr3 = str(eQ - sQ)

								else:
									# compute extra5 bases in case sQ is less
									# than 2 bases bigger than sE
									if (sQ - sE < 2) & (sExtra > 0):
										extra5 = sExtra - sQ + sE

									cds[0] = computeCodons(sQ, mapper, EXON, strand)[1][0]
									cds[1] = computeCodons(eQ, mapper, EXON, strand)[0][1]

									# compute extra3 bases in case eE is less
									# than 2 bases bigger than eQ
									if (eE - eQ < 2) & (eExtra > 0):
										extra3 = eExtra - eE + eQ


						newrow = [cQ, sQ, eQ, hugo, strand, refseq, EXON, intron5, utr5, extra5, cds[0], cds[1], extra3, utr3, intron3, sE, sCpos, eCpos, eE]
						newrow = pd.DataFrame([newrow], columns=cols)

						df = df.append(newrow)


				DF = DF.append(df)

	return(DF)

