#!/usr/bin/python

my=0
ny=0
la=0
si=0
error="error"

for my in range (1,5):
    for ny in range (1, my+1):
        for la in range (1, 5):
            for si in range (1,la):
                print("my,ny,la,si",my,ny,la,si)
                if my>la or (my==la and ny>=si):
                    print (my,ny,la,si)
                elif my<la or (my==la and ny<si):
                    print (la,si,my,ny)
                else: print ('error1')
                if my>=la:
                    if ny>=si:
                        print(my, la, ny, si)
                    elif ny<si:
                        print(my,la,si,ny)
                    else: print ('error2')
                elif my<la:
                    if ny>=si:
                        print(la, my,ny,si)
                    elif ny<si:
                        print(la, my,si,ny)
                    else: print ('error3')
                else: print ('error4')
            si=la
            print("my,ny,la,la",my,ny,la,si)
            if my>la or (my==la and ny>=si):
                print (my,ny,la,si)
            elif my<la or (my==la and ny<si):
                print (la,si,my,ny)
            else: print (error)
            if my>=la:
                if ny>=si:
                    print(my, la, ny, si)
                elif ny<si:
                    print(my,la,si,ny)
                else: print (error)
            elif my<la:
                if ny>=si:
                    print(la, my,ny,si)
                elif ny<si:
                    print(la, my,si,ny)
                else: print (error)
            else: print (error)

