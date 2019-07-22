"""
Read in Surpac code
"""

class Surpac:
    '''Defines a class which reads in the surpac code'''
    def __init__(self, file_name):
        read = open(file_name)
        file = read.readlines()
        read.close()
        header = file[0].strip()
        header_split = header.split(",")
        for i, h in enumerate(header_split):
            header_split[i] = h.strip()
        self.location = header_split[0]
        self.date = header_split[1]
        self.purpose = header_split[2]
        self.memo = header_split[3]
        self.axis = Axis(file[1].strip())
        record = []
        lines = file[2:]
        string_num = None
        temp = None
        self.records = []
        for line in lines:
            string = String(line)
            line_split = line.split(",")
            temp = int(line_split[0])
            if temp != 0:
                if string_num is None or temp == string_num:
                    string_num = temp
                    record.append(string)
                else:
                    raise ValueError('String number changed')
            else:
                if 'END' in string.d[0]:
                    break
                self.records.append(string_record(string_num, record))
                record = []
                string_num = None
                temp = None
    def __str__(self):
        location = str(self.location)
        date = str(self.date)
        purpose = str(self.purpose)
        memo = str(self.memo)
        axis = str(self.axis)
        header = '{}, {}, {}, {}'.format(location, date, purpose, memo)
        record_string = ""
        for i in self.records:
            tog = '{}\n'.format(str(i))
            record_string = record_string + tog
        return '{}\n{}\n{}'.format(header, axis, record_string)

class Axis:
    '''Define a class which reads in the axis line'''
    def __init__(self, axis):
        axis_split = axis.split(",")
        if axis_split[0] != '0':
            raise RuntimeError('First value of axis is not zero')
        self.string = float(axis_split[0])
        self.y_1 = float(axis_split[1])
        self.x_1 = float(axis_split[2])
        self.z_1 = float(axis_split[3])
        self.y_2 = float(axis_split[4])
        self.x_2 = float(axis_split[5])
        self.z_2 = float(axis_split[6])
    def __str__(self):
        string = str(self.string)
        y_1 = str(self.y_1)
        x_1 = str(self.x_1)
        z_1 = str(self.z_1)
        y_2 = str(self.y_2)
        x_2 = str(self.x_2)
        z_2 = str(self.z_2)
        return 'Axis: {}, y1: {}, x1: {}, z1: {}, y2: {},\
             x2: {}, z2: {}'.format(string, y_1, x_1, z_1, y_2, x_2, z_2)

class String:
    '''Class which reads in each line of code and splits it into y, x, z and descriptions'''
    def __init__(self, line):
        line_split = line.split(",")
        self.y = float(line_split[1])
        self.x = float(line_split[2])
        self.z = float(line_split[3])
        self.d = [""]*100
        try:
            for n in range(0, 101):
                self.d[n] = line_split[n+4].strip()
        except IndexError:
            self.length_d = n
    def __str__(self):
        y = str(self.y).strip()
        x = str(self.x).strip()
        z = str(self.z).strip()
        d = list(filter(None, self.d)) 
        description = ''
        for i in d:
            add = '{}, '.format(i)
            description = description + add
        return 'y: {} , x: {} , z: {}, descriptions: {}  '.format(y, x, z, description)

class string_record:
    '''Class which puts all the strings together in one array'''
    def __init__(self, number, strings):
        self.number = number
        self.string_record = strings
    def __str__(self):
        num = str(self.number)
        string_list = ""
        for string in self.string_record:
            line = '{}, {} \n'.format(num, str(string))
            string_list = string_list + line
        return string_list



