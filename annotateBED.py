

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

	"""
	Given a position inside a CDS ("pos") and information about the exon, count
	two pairs of codons:
	- from the start of the CDS to last full codon before "pos" ("UPcds")
	- from the first full codon after "pos" to the end of the CDS ("DOWNcds")
	"""

	# extract exon info from mapper
	sCpos, eCpos = mapper[EXON]['codon_position']
	sC, eC = mapper[EXON]['codon_number']
	sExtra, eExtra = mapper[EXON]['extra_bases']

	### adjust CDS positions and codons numbers in case there are extra bases ###

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


	# get length of CDS segments upstream and downstream of "pos"
	UPl = pos - sCpos
	DOWNl = eCpos - pos

	### count upstream codons, keeping strand into account ###
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


	### count downstream codons, keeping strand into account ###
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



def annotateSymbolTranscriptExonsCodons_02(bed, ucsc):

	"""
	Annotate a BED file (standard, 3 columns) with HUGO symbol, every RefSeq
	transcript ID, every exon number and codons covered. Use a UCSC RefSeq
	reference table.
	"""

	import pandas as pd

	cols = ['chrom', 'start', 'end', 'hugo', 'strand', 'refseq_transcript', 'exon', 'intron5', 'utr5', 'extra_bases5', 'first_complete_codon', 'last_complete_codon', 'extra_bases3', 'utr3', 'intron3', 'exon_start', 'cds_start', 'cds_end', 'exon_end']

	L = []

	# loop over BED intervals (rows)
	for i, row in bed.iterrows():

		if i%10 == 0: print('%s/%s' %(i, len(bed)))

		# get coordinates of QUERY interval
		cQ = row.chrom
		sQ = row.start
		eQ = row.end

		any_hit = False

		# subset REF table for current chromosome
		ucscCHR = ucsc.loc[ucsc.chrom == cQ,].copy()

		# loop over transcripts (rows) that are annotated in this chromosome
		for i2, row2 in ucscCHR.iterrows():

			# get transcript START and END
			sT = row2.txStart
			eT = row2.txEnd

			# if QUERY overlaps with the transcript
			if anyOverlap(sQ, eQ, sT, eT):

				any_hit = True

				# get hugo, transcript ID and strand info
				hugo = row2.name2
				refseq = row2['name']
				strand = row2.strand

				# get every exon START and END position as list of tuples
				starts = row2.exonStarts.split(',')[:-1]
				ends = row2.exonEnds.split(',')[:-1]
				exons = list(zip(starts, ends))

				# get start and end positions of the CDS
				sCDS = int(row2.cdsStart)
				eCDS = int(row2.cdsEnd)

				# create MAPPER for transcript (a dictionary that summarizes exon info)
				mapper = createTranscriptMapper(exons, sCDS, eCDS, strand)

				# track when CDS starts and ends
				cds_status = 'pre_cds'

				# loop over exons:
				for EXON in mapper:

					# get start and end positions of the exon
					sE, eE = mapper[EXON]['position']

					# if QUERY overlaps with the exon
					if anyOverlap(sQ, eQ, sE, eE):

						# get start and end postitions of CDS
						sCpos, eCpos = mapper[EXON]['codon_position']
						# get number of first and last codons
						sC, eC = mapper[EXON]['codon_number']
						# get count of extra bases at start and end of codon
						sExtra, eExtra = mapper[EXON]['extra_bases']

						# adjust codon numbers if extra bases are present
						if sExtra > 0:
							if strand == '+': sC = sC + 1
							else: sC = sC - 1
						if eExtra > 0:
							if strand == '+': eC = eC - 1
							else: eC = eC + 1


						# track when CDS starts and ends
						if (sCpos != 0) & (sE < sCpos):
							cds_status = 'cds'
						if (eCpos != 0) & (eCpos < eE):
							cds_status = 'post_cds'

						# set variables to describe exon
						intron5 = '0'
						utr5 = '0'
						extra5 = '0'
						cds = ['0', '0']
						extra3 = '0'
						utr3 = '0'
						intron3 = '0'

						# if QUERY and the exon overlap like this:
						# EXON			|---------------|
						# QUERY 	|----------------------|
						if (sQ < sE) & (eQ > eE):

							# compute overlap with 5' and 3' introns
							intron5 = str(sE - sQ)
							intron3 = str(eQ - eE)

							# if no CDS is contained in the exon
							if (sCpos, eCpos) == (0, 0):
								# if CDS has not started yet
								if cds_status == 'pre_cds':
									# compute overlap with UTR5
									utr5 = str(eE - sE)
								# if CDS has already ended
								elif cds_status == 'post_cds':
									# compute overlap with UTR3
									utr3 = str(eE - sE)

							# if a CDS is contained in the exon
							else:
								# if the CDS starts in the exon (i.e. the exon contains an UTR5)
								if sE < sCpos:
									# compute overlap with UTR5
									utr5 = str(sCpos - sE)
									# get codons from CDS start to exon end
									cds = (sC, eC)
									# get count of extra bases at 3'
									extra3 = str(eExtra)

								# if the CDS ends in the exon (i.e. the exon contains an UTR3)
								elif eCpos < eE:
									# get count of extra bases at 5'
									extra5 = str(sExtra)
									# get codons from exon start to CDS end
									cds = (sC, eC)
									# compute overlap with UTR3
									utr3 = str(eE - eCpos)

								# if the exon contains only the CDS
								else:
									# get count of extra bases at 5'
									extra5 = str(sExtra)
									# get codons from exon start to exon end
									cds = (sC, eC)
									# get count of extra bases at 3'
									extra3 = str(eExtra)



						# if QUERY and the exon overlap like this:
						# EXON			|------------------|
						# QUERY 	|---------------|
						if (sQ < sE) & (eQ <= eE):

							# compute overlap with 5' intron
							intron5 = str(sE - sQ)

							# if no CDS is contained in the exon
							if (sCpos, eCpos) == (0, 0):
								# if CDS has not started yet
								if cds_status == 'pre_cds':
									# compute overlap with UTR5
									utr5 = str(eQ - sE)
								# if CDS has already ended
								elif cds_status == 'post_cds':
									# compute overlap with UTR3
									utr3 = str(eQ - sE)

							# if a CDS is contained in the exon
							else:
								# if the CDS starts in the exon (i.e. the exon contains an UTR5)
								if sE < sCpos:
									# if there are extra bases 3' and the end of QUERY
									# falls in these extra bases 3' interval
									if (eE - eQ < 2) & (eExtra > 0):
										# adjust count of extra bases at 3'
										extra3 = eExtra - eE + eQ

									# if QUERY does not overlap with CDS
									if eQ <= sCpos:
										# compute overlap with UTR5
										utr5 = str(eQ - sE)
									else:
										# compute overlap with UTR5
										utr5 = str(sCpos - sE)
										# count codons from CDS start to QUERY end
										cds = computeCodons(eQ, mapper, EXON, strand)[0]

								# if the CDS ends in the exon (i.e. the exon contains an UTR3)
								elif eCpos < eE:
									# get count of extra bases at 5'
									extra5 = str(sExtra)
									# if QUERY does not overlap with UTR3
									if eQ <= eCpos:
										# count codons from CDS start to QUERY end
										cds = computeCodons(eQ, mapper, EXON, strand)[0]
									else:
										# get codons from exon start to CDS end
										cds = (sC, eC)
										# compute overlap with UTR3
										utr3 = str(eQ - eCpos)

								# if the exon contains only the CDS
								else:
									# get count of extra bases at 5'
									extra5 = str(sExtra)
									# count codons from exon start to QUERY end
									cds = computeCodons(eQ, mapper, EXON, strand)[0]
									# if there are extra bases 3' and the end of QUERY
									# falls in these extra bases 3' interval
									if (eE - eQ < 2) & (eExtra > 0):
										# adjust count of extra bases at 3'
										extra3 = eExtra - eE + eQ



						# if QUERY and the exon overlap like this:
						# EXON		|------------------|
						# QUERY 	    |------------------|
						if (sQ >= sE) & (eQ > eE):

							# compute overlap with 3' intron
							intron3 = str(eQ - eE)

							# if no CDS is contained in the exon
							if (sCpos, eCpos) == (0, 0):
								# if CDS has not started yet
								if cds_status == 'pre_cds':
									# compute overlap with UTR5
									utr5 = str(eE - sQ)
								# if CDS has already ended
								elif cds_status == 'post_cds':
									# compute overlap with UTR3
									utr3 = str(eE - sQ)

							# if a CDS is contained in the exon
							else:
								# if the CDS starts in the exon (i.e. the exon contains an UTR5)
								if sE < sCpos:
									# get count of extra bases at 3'
									extra3 = str(eExtra)
									# if QUERY overlaps with UTR5
									if sQ < sCpos:
										# compute overlap with UTR5
										utr5 = str(sCpos - sQ)
										# get codons from CDS start to exon end
										cds = (sC, eC)
									# if QUERY does not overlap with UTR5
									else:
										# count codons from QUERY start to exon end
										cds = computeCodons(sQ, mapper, EXON, strand)[1]

								# if the CDS ends in the exon (i.e. the exon contains an UTR3)
								elif eCpos < eE:
									# if there are extra bases 5' and the start of QUERY
									# falls in these extra bases 5' interval
									if (sQ - sE < 2) & (sExtra > 0):
										# adjust count of extra bases at 5'
										extra5 = sExtra - sQ + sE
									# if QUERY overlaps with CDS
									if sQ < eCpos:
										# count codons from QUERY start to CDS end
										cds = computeCodons(sQ, mapper, EXON, strand)[1]
										# compute overlap with UTR3
										utr3 = str(eE - eCpos)
									# if QUERY does not overlap with CDS
									else:
										# compute overlap with UTR3
										utr3 = str(eE - sQ)

								# if the exon contains only the CDS
								else:
									# if there are extra bases 5' and the start of QUERY
									# falls in these extra bases 5' interval
									if (sQ - sE < 2) & (sExtra > 0):
										# adjust count of extra bases at 5'
										extra5 = sExtra - sQ + sE
									# count codons from QUERY start to exon end
									cds = computeCodons(sQ, mapper, EXON, strand)[1]
									# get count of extra bases at 3'
									extra3 = str(eExtra)



						# if QUERY and the exon overlap like this:
						# EXON		|----------------------|
						# QUERY 	    |------------|
						if (sQ >= sE) & (eQ <= eE):

							# if no CDS is contained in the exon
							if (sCpos, eCpos) == (0, 0):
								# if CDS has not started yet
								if cds_status == 'pre_cds':
									# compute overlap with UTR5
									utr5 = str(eQ - sQ)
								# if CDS has already ended
								elif cds_status == 'post_cds':
									# compute overlap with UTR3
									utr3 = str(eQ - sQ)

							# if a CDS is contained in the exon
							else:
								# if the CDS starts in the exon (i.e. the exon contains an UTR5)
								if sE < sCpos:
									# if QUERY only overlaps with UTR5
									if eQ < sCpos:
										# compute overlap with UTR5
										utr5 = str(eQ - sQ)
									# if QUERY overlaps with both UTR5 and CDS
									elif sQ < sCpos < eQ:
										# compute overlap with UTR5
										utr5 = str(sCpos - sQ)
										# count codons from CDS start to QUERY end
										cds = computeCodons(eQ, mapper, EXON, strand)[0]
										# if there are extra bases 3' and the end of QUERY
										# falls in these extra bases 3' interval
										if (eE - eQ < 2) & (eExtra > 0):
											# adjust count of extra bases at 3'
											extra3 = eExtra - eE + eQ
									# if QUERY only overlaps with CDS
									elif sCpos < sQ:
										# count codon at QUERY start
										cds[0] = computeCodons(sQ, mapper, EXON, strand)[1][0]
										# count codon at QUERY end
										cds[1] = computeCodons(eQ, mapper, EXON, strand)[0][1]
										# if there are extra bases 3' and the end of QUERY
										# falls in these extra bases 3' interval
										if (eE - eQ < 2) & (eExtra > 0):
											# adjust count of extra bases at 3'
											extra3 = eExtra - eE + eQ

								# if the CDS ends in the exon (i.e. the exon contains an UTR3)
								elif eCpos < eE:
									# if QUERY overlaps with CDS only
									if eQ < eCpos:
										# count codon at QUERY start
										cds[0] = computeCodons(sQ, mapper, EXON, strand)[1][0]
										# count codon at QUERY end
										cds[1] = computeCodons(eQ, mapper, EXON, strand)[0][1]
										# if there are extra bases 5' and the start of QUERY
										# falls in these extra bases 5' interval
										if (sQ - sE < 2) & (sExtra > 0):
											# adjust count of extra bases at 5'
											extra5 = sExtra - sQ + sE
									# if QUERY overlaps with both CDS and UTR3
									elif sQ < eCpos < eQ:
										# if there are extra bases 5' and the start of QUERY
										# falls in these extra bases 5' interval
										if (sQ - sE < 2) & (sExtra > 0):
											# adjust count of extra bases at 5'
											extra5 = sExtra - sQ + sE
										# count codons from QUERY start to CDS end
										cds = computeCodons(sQ, mapper, EXON, strand)[1]
										# compute overlap with UTR3
										utr3 = str(eQ - eCpos)
									# if QUERY overlaps with UTR3 only
									elif eCpos < sQ:
										# compute overlap with UTR3
										utr3 = str(eQ - sQ)

								# if the exon contains only the CDS
								else:
									# if there are extra bases 5' and the start of QUERY
									# falls in these extra bases 5' interval
									if (sQ - sE < 2) & (sExtra > 0):
										# adjust count of extra bases at 5'
										extra5 = sExtra - sQ + sE
									# count codon at QUERY start
									cds[0] = computeCodons(sQ, mapper, EXON, strand)[1][0]
									# count codon at QUERY end
									cds[1] = computeCodons(eQ, mapper, EXON, strand)[0][1]
									# if there are extra bases 3' and the end of QUERY
									# falls in these extra bases 3' interval
									if (eE - eQ < 2) & (eExtra > 0):
										# adjust count of extra bases at 3'
										extra3 = eExtra - eE + eQ

						# create new row for final DF and append to growing list
						newrow = [cQ, sQ, eQ, hugo, strand, refseq, EXON, intron5, utr5, extra5, cds[0], cds[1], extra3, utr3, intron3, sE, sCpos, eCpos, eE]
						# if strand is negative, switch the order between 5' and 3' elements
						if strand == '-':
							newrow = [cQ, sQ, eQ, hugo, strand, refseq, EXON, intron3, utr3, extra3, cds[1], cds[0], extra5, utr5, intron5, sE, sCpos, eCpos, eE]
						L.append(newrow)

		if not any_hit:
			newrow = [cQ, sQ, eQ] + ['' for i in range(len(cols[3:]))]
			L.append(newrow)

	# transform to DF
	DF = pd.DataFrame(L, columns=cols)
	return(DF)
