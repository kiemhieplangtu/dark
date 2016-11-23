print 'hi!'
month = int(raw_input('Please enter month: '))

print 'Month: ', month

if month<0 or month>12:
	print 'Try again!'
	month = int(raw_input('Please enter month: '))
