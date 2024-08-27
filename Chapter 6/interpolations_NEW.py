import numpy as np
import warnings

def quadratic_function(z, psi2, psi1, psi4):
    a = 2*(psi1 + psi2 -2*psi4)
    b =  (4*psi4 - 3*psi2 -psi1)
    c = psi2
    return a * z**2 + b * z + c

def triangle_area(x_coords, y_coords):
    # Extract coordinates of the vertices
    x1= x_coords[0]
    x2= x_coords[1]
    x3= x_coords[2]
    y1= y_coords[0]
    y2= y_coords[1]
    y3= y_coords[2]

    # Calculate the area using the formula: 0.5 * |(x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2))|
    area = 0.5 * abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))) 
    return area

warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in sqrt')

def IntersectionLength(psi2,psi1,psi4,positive=True):# (begin,end,mid) #returns the ratio between the cut part and the side
    a = 2*(psi1 + psi2 -2*psi4)
    b =  (4*psi4 - 3*psi2 -psi1)
    c = psi2
    discriminant = b**2 - 4 * a * c
    if (discriminant) < 0 or (a == 0):
        return 1
    if positive:
        L_bar_ratio =( -b+np.sqrt(b**2-4*a*c) ) / (2*a)  
    else:
        L_bar_ratio =( -b-np.sqrt(b**2-4*a*c) ) / (2*a)
    if not np.isfinite(L_bar_ratio) or L_bar_ratio == 0:
        return 1
    return L_bar_ratio

def GetVolumeFraction2(psi0_5,x_coords,y_coords,EPS):
    s = 0
    #Case 0  
    print("psi0_5",psi0_5)
    psi0_5 = np.array(psi0_5)

    if all(p < -EPS for p in psi0_5):
        return 1
    elif all(p > EPS for p in psi0_5):
        return 0

    # #Constructing the linear system AX=B
    # A = np.zeros((6, 6))

    #         # Loop over each point

    # for k in range(6):
    #     x = x_coords[k]
    #     y = y_coords[k]
                
    #         # Construct the row for A
    #     A[k] = [x**2, y**2, x*y, x, y, 1]

   ## Cases determinators
    positive_count = np.sum(psi0_5 > EPS)
    negative_count = np.sum(psi0_5 < -EPS)
    zero_count = np.sum(abs(psi0_5) < EPS)
    print("[#+ , #-, #0]",positive_count,negative_count,zero_count)

    #Case1 , all non zero , one has opposite sign
    if (positive_count == 5 and negative_count == 1) or (positive_count == 1 and negative_count == 5):
        negative_index =  np.where(psi0_5 < -EPS)[0] if (positive_count == 5 and negative_count == 1) else np.where(psi0_5  > EPS)[0]            
        
        if negative_index[0] in [2, 0, 1]: ## Subcase 1.1 solo negative on vertex 
            line_pairs = [(2, 1, 4, 0, 3), (0, 2, 3, 1, 5), (1, 0, 5, 2, 4)]
            for i, pair in enumerate(line_pairs):
                if negative_index[0] == pair[0]:
                    psi_start, psi_end1, psi_mid1 = psi0_5[pair[0]], psi0_5[pair[1]], psi0_5[pair[2]]
                    psi_end2, psi_mid2 = psi0_5[pair[3]], psi0_5[pair[4]]

                    l_11 = IntersectionLength(psi_start, psi_end1, psi_mid1)
                    l_12 = IntersectionLength(psi_start, psi_end1, psi_mid1, False)
                    l_1 = l_11 if 0 < l_11 < 0.5 else l_12
                    
                    l_21 = IntersectionLength(psi_start, psi_end2, psi_mid2)
                    l_22 = IntersectionLength(psi_start, psi_end2, psi_mid2, False)
                    l_2 = l_21 if 0 < l_21 < 0.5 else l_22           
                    s = l_1 * l_2
                    break # Changed to return result immediately
                    
        elif negative_index[0] in [4, 3, 5]: ## Subcase 1.2 solo negative on midpoint 
            side_pairs = [(2, 0), (2, 1), (1, 0)]
            for i, pair in enumerate(side_pairs):
                if negative_index[0] == (i + 3):
                    L_side = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)
                    psi_start = psi0_5[pair[0]]
                    psi_end = psi0_5[pair[-1]]
                    psi_mid = psi0_5[negative_index[0]]

                    l_11 = IntersectionLength(psi_start, psi_end, psi_mid)
                    l_12 = IntersectionLength(psi_start, psi_end, psi_mid, False)
                    l_1 = l_11 * L_side
                    l_2 = l_12 * L_side

                    area = triangle_area(x_coords, y_coords)
                    z_values = np.linspace(l_1, l_2, num=10)
                    function_values = quadratic_function(z_values, psi_start, psi_end, psi_mid)
                    area_under = np.abs(np.trapz(function_values, x=z_values))
                    s = area_under / area
                    break
                    
        if (positive_count == 1 and negative_count == 5):
            s = 1-s
        return s
        
    #Case2, all non zeros , two has opposite  sign 
    elif (positive_count == 4 and negative_count == 2) or (positive_count == 2 and negative_count == 4): 
        #Ver and consecative () ,ver and opposite mid ,2 mids , 2 vers 
                                                    ## 2 with 3 or 4,, zero with 3 or 5 ,, one with 4 or 5  
        negative_index =  np.where(psi0_5 < -EPS)[0] if (positive_count == 4 and negative_count == 2) else np.where(psi0_5 > EPS)[0]

        # Defining pair sets for all possible combinations (4 Subcases) #ASK WHICH ARE IMPOSSIBLE 
    
        vertex_midpoint_con_pairs = [
                (2, 3), (2, 4),  
                (0, 3), (0, 5),  
                (1, 4), (1, 5)   
            ]
        vertex_midpoint_oppo_pairs = [
                (2, 5),  
                (0, 4),  
                (1, 3)
        ]
        vertex_vertex_pairs = [
                (2, 1), (2, 0),  
                (1, 0)  
            ]
        midpoint_midpoint_pairs = [
                (3, 4), (3, 5),  
                (4, 5)  
            ]
    
        for pair in vertex_midpoint_con_pairs:
            if set(negative_index) == set(pair):  # Check if the set of negative indices matches the current pair
                vertex, midpoint = pair  # Unpack the vertex and midpoint indices from the pair
                psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values for vertices 0, 1, and 2
                psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5]  # Extract level-set values for midpoints 3, 4, and 5
               
                if vertex == 2 and midpoint == 3:  # If vertex 2 is negative
                    
                    l_11 = IntersectionLength(psi_v2, psi_v1, psi_m4) 
                    l_12 = IntersectionLength(psi_v2, psi_v1, psi_m4, False) 
                    l_1 = l_11 if l_11 < 0.5 else l_12  
    
                    l_21 = IntersectionLength(psi_v2, psi_v0, psi_m3)  
                    l_22 = IntersectionLength(psi_v2, psi_v0, psi_m3, False) 
                    l_2 = l_21 if 0.5 < l_21  else l_22  
                elif vertex ==2 and midpoint ==4:
                    l_11 = IntersectionLength(psi_v2, psi_v1, psi_m4)  
                    l_12 = IntersectionLength(psi_v2, psi_v1, psi_m4, False)   
                    l_1 = l_11 if l_11 > 0.5 else l_12   
    
                    l_21 = IntersectionLength(psi_v2, psi_v0, psi_m3)   
                    l_22 = IntersectionLength(psi_v2, psi_v0, psi_m3, False)   
                    l_2 = l_21 if 0 < l_21 < 0.5 else l_22  
                elif vertex ==0 and midpoint ==3:
                                    
                    l_11 = IntersectionLength(psi_v0, psi_v1, psi_m5)  
                    l_12 = IntersectionLength(psi_v0, psi_v1, psi_m5, False)   
                    l_1 = l_11 if l_11 > 0.5 else l_12   
    
                    l_21 = IntersectionLength(psi_v0, psi_v2, psi_m3)  
                    l_22 = IntersectionLength(psi_v0, psi_v2, psi_m3, False)  
                    l_2 = l_21 if 0 < l_21 < 0.5 else l_22  
                elif vertex ==0 and midpoint ==5:
                                                    
                    l_11 = IntersectionLength(psi_v0, psi_v1, psi_m5)  
                    l_12 = IntersectionLength(psi_v0, psi_v1, psi_m5, False) 
                    l_1 = l_11 if 0 < l_11 < 0.5 else l_12  #
    
                    l_21 = IntersectionLength(psi_v0, psi_v2, psi_m3) 
                    l_22 = IntersectionLength(psi_v0, psi_v2, psi_m3, False)  
                    l_2 = l_21 if 0.5 > l_21  else l_22 
                elif vertex ==1 and midpoint ==3:
                                                                    
                    l_11 = IntersectionLength(psi_v0, psi_v1, psi_m5)  
                    l_12 = IntersectionLength(psi_v0, psi_v1, psi_m5, False) 
                    l_1 = l_11 if 0 < l_11 < 0.5 else l_12  #
    
                    l_21 = IntersectionLength(psi_v0, psi_v2, psi_m3) 
                    l_22 = IntersectionLength(psi_v0, psi_v2, psi_m3, False)  
                    l_2 = l_21 if 0.5 > l_21  else l_22 
                elif vertex ==1 and midpoint ==4:
                    
                    l_11 = IntersectionLength(psi_v1, psi_v0, psi_m5)  
                    l_12 = IntersectionLength(psi_v1, psi_v0, psi_m5, False) 
                    l_1 = l_11 if 0 < l_11 < 0.5 else l_12  #
    
                    l_21 = IntersectionLength(psi_v1, psi_v2, psi_m4) 
                    l_22 = IntersectionLength(psi_v1, psi_v2, psi_m4, False)  
                    l_2 = l_21 if 0.5 > l_21  else l_22 
                else:
                    break
                s = l_1*l_2
        for pair in vertex_midpoint_oppo_pairs:
            if set(negative_index) == set(pair):  # Check if the set of negative indices matches the current pair
                vertex, midpoint = pair  # Unpack the vertex and midpoint indices from the pair
                psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values for vertices 0, 1, and 2
                psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5]  # Extract level-set values for midpoints 3, 4, and 5
                if vertex == 2: # If vertex 2 is negative >>> Mid 5 is negtive
                    
                    l_11 = IntersectionLength(psi_v2, psi_v1, psi_m4) 
                    l_12 = IntersectionLength(psi_v2, psi_v1, psi_m4, False) 
                    l_1 = l_11 if l_11 < 0.5 else l_12  
                    l_21 = IntersectionLength(psi_v2, psi_v0, psi_m3) 
                    l_22 = IntersectionLength(psi_v2, psi_v0, psi_m3, False) 
                    l_2 = l_21 if l_21 < 0.5 else l_22  
                    area = triangle_area(x_coords, y_coords)
                    area_vertex = l_1*l_2*area

                    l_31 = IntersectionLength(psi_v0, psi_v1, psi_m5) 
                    l_32 = IntersectionLength(psi_v0, psi_v1, psi_m5, False) 
                    L_side = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)

                    l_31 = l_31*L_side
                    l_32 = l_32*L_side

                
                    z_values = np.linspace(l_31, l_32, num=10)
                    function_values = quadratic_function(z_values, psi_v0, psi_v1, psi_m5)
                    area_under = np.abs(np.trapz(function_values, x=z_values))
                    s = (area_under+area_vertex)/(area)
                
                elif vertex == 1: # If vertex 1 is negative >>> Mid 3 is negtive
                    
                    l_11 = IntersectionLength(psi_v1, psi_v2, psi_m4) 
                    l_12 = IntersectionLength(psi_v1, psi_v2, psi_m4, False) 
                    l_1 = l_11 if l_11 < 0.5 else l_12  
                    l_21 = IntersectionLength(psi_v1, psi_v0, psi_m5) 
                    l_22 = IntersectionLength(psi_v1, psi_v0, psi_m5, False) 
                    l_2 = l_21 if l_21 < 0.5 else l_22  
                    area = triangle_area(x_coords, y_coords)
                    area_vertex = l_1*l_2*area

                    l_31 = IntersectionLength(psi_v0, psi_v2, psi_m3) 
                    l_32 = IntersectionLength(psi_v0, psi_v2, psi_m3, False) 

                    L_side = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)

                    l_31 = l_31*L_side
                    l_32 = l_32*L_side

                    z_values = np.linspace(l_31, l_32, num=10)
                    
                    function_values = quadratic_function(z_values, psi_v0, psi_v2, psi_m3)
                    area_under = np.abs(np.trapz(function_values, x=z_values))
                    s = (area_under+area_vertex)/(area)
                elif vertex == 0: # If vertex 0 is negative >>> Mid 4 is negtive
                    
                    l_11 = IntersectionLength(psi_v0, psi_v2, psi_m3) 
                    l_12 = IntersectionLength(psi_v0, psi_v2, psi_m3, False) 
                    l_1 = l_11 if l_11 < 0.5 else l_12  
                    l_21 = IntersectionLength(psi_v0, psi_v1, psi_m5) 
                    l_22 = IntersectionLength(psi_v0, psi_v1, psi_m5, False) 
                    l_2 = l_21 if l_21 < 0.5 else l_22  
                    area = triangle_area(x_coords, y_coords)
                    area_vertex = l_1*l_2*area

                    l_31 = IntersectionLength(psi_v1, psi_v2, psi_m4) 
                    l_32 = IntersectionLength(psi_v1, psi_v2, psi_m4, False) 


                    L_side = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)

                    l_31 = l_31*L_side
                    l_32 = l_32*L_side

                
                
                    z_values = np.linspace(l_31, l_32, num=10)
                    
                    function_values = quadratic_function(z_values, psi_v1, psi_v2, psi_m3)
                    area_under = np.abs(np.trapz(function_values, x=z_values))
                    s = (area_under+area_vertex)/(area)
                    
        for pair in vertex_vertex_pairs: #(0,1)(0,2)(1,2)
            if set(negative_index) == set(pair):  # Check if the set of negative indices matches the current pair 
                #if vertex vertex, then one mid point is shared and one is not 
                
                psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values for vertices 0, 1, and 2
                psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5]  # Extract level-set values for midpoints 3, 4, and 5

                if set(negative_index) == set((1,0)):
                    psi_v1st = psi_v0
                    psi_v2nd = psi_v1
                    psi_m_shared = psi_m5
                    psi_m_v1st = psi_m3
                    psi_m_v2nd = psi_m4
                    psi_v_out = psi_v2
                elif  set(negative_index) == set((2,0)):
                    psi_v1st = psi_v0
                    psi_v2nd = psi_v2
                    psi_m_shared = psi_m3
                    psi_m_v1st = psi_m5
                    psi_m_v2nd = psi_m4
                    psi_v_out = psi_v1
                elif set(negative_index) == set((1,2)): 
                    psi_v1st = psi_v1
                    psi_v2nd = psi_v2
                    psi_m_shared = psi_m4
                    psi_m_v1st = psi_m5
                    psi_m_v2nd = psi_m3
                    psi_v_out = psi_v0


                #we need four lengths so 8 calls
                l_11 = IntersectionLength(psi_v1st, psi_v2nd, psi_m_shared)
                l_12 = IntersectionLength(psi_v1st, psi_v2nd, psi_m_shared,False)
                l_1 = l_11 if (l_11>0 and l_11 < 0.5) else l_12 

                l_21 = IntersectionLength(psi_v2nd, psi_v1st, psi_m_shared)
                l_22 = IntersectionLength(psi_v2nd, psi_v1st, psi_m_shared,False)
                l_2 = l_21 if (l_21>0 and l_21 < 0.5) else l_22 

                l_31 = IntersectionLength(psi_v1st, psi_v_out, psi_m_v1st)
                l_32 = IntersectionLength(psi_v1st, psi_v_out, psi_m_v1st,False)
                l_3 = l_31 if (l_31>0 and l_31 < 0.5) else l_32 


                l_41 = IntersectionLength(psi_v2nd, psi_v_out, psi_m_v2nd)
                l_42 = IntersectionLength(psi_v2nd, psi_v_out, psi_m_v2nd,False)
                l_4 = l_41 if (l_41>0 and l_41 < 0.5) else l_42 

                s = (l_1*l_3)+(l_2*l_4)
                
                # l_11 = IntersectionLength(psi_v0, psi_v2, psi_m3) 
                # l_12 = IntersectionLength(psi_v0, psi_v2, psi_m3, False) 
                # l_1 = l_11 if l_11 < 0.5 else l_12  
                # l_21 = IntersectionLength(psi_v0, psi_v1, psi_m5) 
                # l_22 = IntersectionLength(psi_v0, psi_v1, psi_m5, False) 
                # l_2 = l_21 if l_21 < 0.5 else l_22  
                # s0 = l_1*l_2

                # l_11 = IntersectionLength(psi_v1, psi_v0, psi_m5) 
                # l_12 = IntersectionLength(psi_v1, psi_v0, psi_m5, False) 
                # l_1 = l_11 if l_11 < 0.5 else l_12  
                # l_21 = IntersectionLength(psi_v1, psi_v2, psi_m4) 
                # l_22 = IntersectionLength(psi_v1, psi_v2, psi_m4, False) 
                # l_2 = l_21 if l_21 < 0.5 else l_22  
                # s1 = l_1*l_2

                # l_11 = IntersectionLength(psi_v2, psi_v0, psi_m3) 
                # l_12 = IntersectionLength(psi_v2, psi_v0, psi_m3, False) 
                # l_1 = l_11 if l_11 < 0.5 else l_12  
                # l_21 = IntersectionLength(psi_v2, psi_v1, psi_m4) 
                # l_22 = IntersectionLength(psi_v2, psi_v1, psi_m4, False) 
                # l_2 = l_21 if l_21 < 0.5 else l_22  
                # s2 = l_1*l_2
                
                # if vertex1 == 0 and vertex2 == 1:
                #     s = s0 + s1
                # elif vertex1 == 0 and vertex2 == 2:
                #     s = s0 + s2
                # elif vertex1 == 1 and vertex2 == 2:
                #     s = s1 + s2
        for pair in midpoint_midpoint_pairs: #(3,4)(3,5)(4,5)
            if set(negative_index) == set(pair):  
                midpoint1, midpoint2 = pair  # Unpack the vertex and midpoint indices from the pair
                psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values for vertices 0, 1, and 2
                psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5]  # Extract level-set values for midpoints 3, 4, and 5
               #Calculating the area and area under for each 
                L_01 = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)
                L_02 = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)
                L_12 = np.sqrt((x_coords[pair[0]] - x_coords[pair[1]])**2 + (y_coords[pair[0]] - y_coords[pair[1]])**2)
                
                #for point 3
                l_31 = IntersectionLength(psi_v0, psi_v1, psi_m5) * L_01
                l_32 = IntersectionLength(psi_v0, psi_v1, psi_m5, False) *L_01 
                z_values = np.linspace(l_31, l_32, num=10)
                function_values = quadratic_function(z_values, psi_v0, psi_v1, psi_m5)
                area3= np.abs(np.trapz(function_values, x=z_values))

                #For point 4
                l_41 = IntersectionLength(psi_v1, psi_v2, psi_m4) * L_12
                l_42 = IntersectionLength(psi_v1, psi_v2, psi_m4, False) *L_12 
                z_values = np.linspace(l_41, l_42, num=10)
                function_values = quadratic_function(z_values, psi_v1, psi_v2, psi_m4)
                area4= np.abs(np.trapz(function_values, x=z_values))

                #For point 5 
                l_51 = IntersectionLength(psi_v0, psi_v2, psi_m3) * L_02
                l_52 = IntersectionLength(psi_v0, psi_v2, psi_m3, False) *L_02 
                z_values = np.linspace(l_51, l_52, num=10)
                function_values = quadratic_function(z_values, psi_v0, psi_v2, psi_m3)
                area5= np.abs(np.trapz(function_values, x=z_values))

                area = triangle_area(x_coords, y_coords)

                if midpoint1 == 3 and midpoint2 ==4:
                    s = (area3+area4)/area
                elif midpoint1 == 3 and midpoint2 ==5:
                    s = (area3+area5)/area
                elif midpoint1 == 4 and midpoint2 ==5:
                    s = (area4+area5)/area

        if (positive_count == 2 and negative_count == 4):
            s = 1-s
        return s    
    #Case 3, all non zeros, three and three                     
    elif positive_count == 3 and negative_count == 3: # Corner or line , total 6 case options
                                                    ## Corners : 2-3-4 __ 0-3-5  __1-5-4
                                                    ## Sides : 241 __ 150__230 >>>> side5 = 1-corner10
                                                                                    # side6 = 1-corner9
                                                                                    # side7 = 1- corner8 
        negative_index =  np.where(psi0_5 < -EPS)[0]



        psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values for vertices 0, 1, and 2
        psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5]  # Extract level-set values for midpoints 3, 4, and 5
        
        sum_negative = np.sum(negative_index) #this array has the sum of index numbers , the number means nothing but it relates the cases 
        if sum_negative == 5 or sum_negative == 10: #side_230(5)=1-corner154(10)
            l_11 = IntersectionLength(psi_v1,psi_v2,psi_m4)
            l_12 = IntersectionLength(psi_v1,psi_v2,psi_m4,False)
            l_1  = l_11 if l_11>0.5 else l_12

            l_21 = IntersectionLength(psi_v1,psi_v0,psi_m5)
            l_22 = IntersectionLength(psi_v1,psi_v0,psi_m5,False)

            l_2 = l_21 if l_21 > 0.5 else l_22
      
            s = l_1*l_2
            if sum_negative == 5:
                s = 1-s 
        elif sum_negative == 6 or sum_negative == 9: #side_150 #corner243
            l_11 = IntersectionLength(psi_v2,psi_v0,psi_m3)
            l_12 = IntersectionLength(psi_v2,psi_v0,psi_m3,False)
            l_1  = l_11 if l_11>0.5 else l_12

            l_21 = IntersectionLength(psi_v2,psi_v1,psi_m4)
            l_22 = IntersectionLength(psi_v2,psi_v1,psi_m4,False)

            l_2 = l_21 if l_21 > 0.5 else l_22
      
            s  = l_1*l_2
            if sum_negative == 6:
                s = 1-s 
        elif sum_negative == 7 or sum_negative == 8: #side_241 #corner 3o5
            l_11 = IntersectionLength(psi_v0,psi_v2,psi_m3)
            l_12 = IntersectionLength(psi_v0,psi_v2,psi_m3,False)
            l_1  = l_11 if l_11>0.5 else l_12

            l_21 = IntersectionLength(psi_v0,psi_v1,psi_m5)
            l_22 = IntersectionLength(psi_v0,psi_v1,psi_m5,False)

            l_2 = l_21 if l_21 > 0.5 else l_22
      
            s  = l_1*l_2
            if sum_negative == 7:
                s = 1-s 
        return s 
    #Case 4 
    elif zero_count == 1 and negative_count == 0:
        return 0
    elif zero_count == 1 and positive_count == 0:
        return 1

    #Case 5, one zero, one opposite 
    elif zero_count == 1 and (negative_count == 1 or positive_count == 1 ): #if the zero is on a vertex , the negtive has to be on a concecative midpoint (6cases)
                                                                                            # V2 with M5 or M4 _ V0 with M5 or M3 _ V1 with M4 or M3  
                                                  #if the zero is on a midpoint , the negative hast to be on one of the side's vertices (6cases) 
                                                                                            # M3 woth V2 or V0 , M5 with V1 or V0 , M4 with V2 or V1
        
        negative_index =  np.where(psi0_5 < -EPS)[0] if negative_count ==1 else np.where(psi0_5 < -EPS)[0] 


        zero_index = np.where(abs(psi0_5) < EPS)[0]
            
        vertex_midpoint_con_pairs = [
                (2, 3), (2, 4),  
                (0, 3), (0, 5),  
                (1, 4), (1, 5)   ]
        midpoint_vertex_con_pairs = [
                (5, 0), (5, 1),  
                (4, 1), (4, 2),  
                (3, 0), (3, 2)      ]
        
        psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values 
        psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5] 
         # SUBCASE 5.2
        if (zero_index[0] == 5 and negative_index[0] == 0   ): #the or is probably redundant 
            l_11 = IntersectionLength(psi_v0,psi_v1,psi_m5)
            l_12 = IntersectionLength(psi_v0,psi_v1,psi_m5,False) 
            l_1 = l_11 if l_11 <1/2 else l_12
            s = l_1/2
        elif (zero_index[0] == 5 and negative_index[0] == 2  ):
            l_11 = IntersectionLength(psi_v2,psi_v1,psi_m4)
            l_12 = IntersectionLength(psi_v2,psi_v1,psi_m4,False) 
            l_1 = l_11 if l_11 <1/2 else l_12
            s = l_1/2
        elif (zero_index[0] == 4 and negative_index[0] == 1  ):
            l_11 = IntersectionLength(psi_v1,psi_v0,psi_m5)
            l_12 = IntersectionLength(psi_v1,psi_v0,psi_m5,False) 
            l_1 = l_11 if l_11 <1/2 else l_12
            s = l_1/2
        elif (zero_index[0] == 4 and negative_index[0] == 2  ):
            l_11 = IntersectionLength(psi_v2,psi_v0,psi_m3)
            l_12 = IntersectionLength(psi_v2,psi_v0,psi_m3,False) 
            l_1 = l_11 if l_11 <1/2 else l_12
            s = l_1/2
        elif (zero_index[0] == 3 and negative_index[0] == 0  ):
            l_11 = IntersectionLength(psi_v2,psi_v1,psi_m4)
            l_12 = IntersectionLength(psi_v2,psi_v1,psi_m4,False) 
            l_1 = l_11 if l_11 <1/2 else l_12
            s = l_1/2
        elif (zero_index[0] == 3 and negative_index[0] == 1  ):
            l_11 = IntersectionLength(psi_v2,psi_v1,psi_m4)
            l_12 = IntersectionLength(psi_v2,psi_v1,psi_m4,False) 
            l_1 = l_11 if l_11 <1/2 else l_12
            s = l_1/2
        # SUBCASE 5.1
        elif (zero_index[0] == 2 and negative_index[0] == 3)   :
            L_side = np.sqrt((x_coords[2]**2-x_coords[0]**2)+(y_coords[2]**2-y_coords[0]**2)) #length of side 2-0
            l_11 = IntersectionLength(psi_v2,psi_v0,psi_m3)
            l_12 = IntersectionLength(psi_v2,psi_v0,psi_m3,False)
            
            l_1 = l_11 if l_11 >1/2 else l_12
            l_1 = l_1*L_side 
            
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, l_1, num=10)
            function_values = quadratic_function(z_values, psi_v2, psi_v0, psi_m5)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area
            
        elif (zero_index[0] == 2 and (negative_index[0] == 4 ) ):
            L_side = np.sqrt((x_coords[2]**2-x_coords[1]**2)+(y_coords[2]**2-y_coords[1]**2)) #length of side 2-1
            l_11 = IntersectionLength(psi_v0,psi_v1,psi_m4)
            l_12 = IntersectionLength(psi_v0,psi_v1,psi_m4,False)
            l_1 = l_11 if l_11 >1/2 else l_12
            l_1 = l_1*L_side 
            


            l_1 = l_1*L_side 
            
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, l_1, num=10)
            function_values = quadratic_function(z_values, psi_v2, psi_v1, psi_m4)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area
            
        elif (zero_index[0] == 0 and negative_index[0] == 3):
            L_side = np.sqrt((x_coords[0]**2-x_coords[2]**2)+(y_coords[0]**2-y_coords[2]**2)) #length of side 0-2
            l_11 = IntersectionLength(psi_v0,psi_v2,psi_m3)
            l_12 = IntersectionLength(psi_v0,psi_v2,psi_m3,False)
            l_1 = l_11 if l_11 >1/2 else l_12
            l_1 = l_1*L_side 
            


            l_1 = l_1*L_side 
            
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, l_1, num=10)
            function_values = quadratic_function(z_values, psi_v0, psi_v1, psi_m3)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area    
        
        elif (zero_index[0] == 0 and negative_index[0] == 5):
            L_side = np.sqrt((x_coords[0]**2-x_coords[1]**2)+(y_coords[0]**2-y_coords[1]**2)) #length of side 0-1
            l_11 = IntersectionLength(psi_v0,psi_v1,psi_m5)
            l_12 = IntersectionLength(psi_v0,psi_v1,psi_m5,False)
            l_1 = l_11 if l_11 >1/2 else l_12
            l_1 = l_1*L_side 
            


            l_1 = l_1*L_side 
            
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, l_1, num=10)
            function_values = quadratic_function(z_values, psi_v0, psi_v1, psi_m5)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area      
    
        elif (zero_index[0] == 1 and negative_index[0] == 4):
            L_side = np.sqrt((x_coords[1]**2-x_coords[2]**2)+(y_coords[1]**2-y_coords[2]**2)) #length of side 1-2
            l_11 = IntersectionLength(psi_v1,psi_v2,psi_m4)
            l_12 = IntersectionLength(psi_v1,psi_v2,psi_m4,False)
            l_1 = l_11 if l_11 >1/2 else l_12
            l_1 = l_1*L_side 
            
            
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, l_1, num=10)
            function_values = quadratic_function(z_values, psi_v1, psi_v2, psi_m4)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area

        elif (zero_index[0] == 1 and negative_index[0] == 3):
            L_side = np.sqrt((x_coords[1]**2-x_coords[0]**2)+(y_coords[1]**2-y_coords[0]**2)) #length of side 1-0
            l_11 = IntersectionLength(psi_v1,psi_v0,psi_m3)
            l_12 = IntersectionLength(psi_v1,psi_v0,psi_m3,False)
            l_1 = l_11 if l_11 >1/2 else l_12
            l_1 = l_1*L_side 
            


            l_1 = l_1*L_side 
            
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, l_1, num=10)
            function_values = quadratic_function(z_values, psi_v1, psi_v2, psi_m4)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area
        if positive_count == 1:
            s = 1-s
        return s 
        #Case 6 , one zero two opposites cant happen ?DOCH  > (#+ , #-, #0)= (2 3 1), ii,el,s= 3630 (V1878, V1911, V1905) None
    elif zero_count == 1 and (negative_count == 2 or positive_count == 2 ):
        return 0 
        #Case 6' , two zeros one opposite > Cant have consecative zeros , either two vertices or two midpoints, each two zeros allow a negative only between them 
    elif zero_count == 2 and (negative_count == 1 or positive_count == 1 ): 
        
        negative_index =  np.where(psi0_5 < -EPS)[0] if negative_count ==1 else np.where(psi0_5 < -EPS)[0] 

        zero_index = np.where(abs(psi0_5) < EPS)[0]
            
        vertex_pairs = [ #Subcase6.1  >>> area under lune
                (0, 1, 5), (0, 2, 3),  
                (1, 2, 4) 
        ]
        midpoint_pairs = [#Subcase6.2 >>>> corner area
            (3, 4, 2),  (3, 5, 0), 
            (4, 5, 1)  
        ]

        
        psi_v0, psi_v1, psi_v2 = psi0_5[0], psi0_5[1], psi0_5[2]  # Extract level-set values 
        psi_m3, psi_m4, psi_m5 = psi0_5[3], psi0_5[4], psi0_5[5] 
       #SUBCASE 6.1
        if negative_index == 5:
            L_side = np.sqrt((x_coords[0]**2-x_coords[1]**2)+(y_coords[0]**2-y_coords[1]**2)) #length of side 0-1
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, L_side, num=10)
            function_values = quadratic_function(z_values, 0, 0, psi_m5)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area
        elif negative_index == 4:
            L_side = np.sqrt((x_coords[1]**2-x_coords[2]**2)+(y_coords[1]**2-y_coords[2]**2)) #length of side 1-2
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, L_side, num=10)
            function_values = quadratic_function(z_values, 0, 0, psi_m4)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area
        elif negative_index == 3:
            L_side = np.sqrt((x_coords[0]**2-x_coords[2]**2)+(y_coords[0]**2-y_coords[2]**2)) #length of side 0-2
            area = triangle_area(x_coords, y_coords)
            z_values = np.linspace(0, L_side, num=10)
            function_values = quadratic_function(z_values, 0, 0, psi_m4)
            area_under = np.abs(np.trapz(function_values, x=z_values))
            s = area_under/area
        #SUBCASE6.2
        elif (negative_index ==2) or (negative_index ==1) or (negative_index ==0) :
            s = 1/4

        if positive_count == 1:
            s = 1-s
        return s

from ngsolve import *

def p2_InterpolateLevelSetToElems(levelset_p2, val1, val2, func_p0, mesh, EPS):
    s  = 0
    nCutElems = 0
    fes= H1(mesh, order=2)
    #for el in mesh.Elements():  
    for ii,el in enumerate(mesh.Elements()): # in  each of (mesh.ne) elements , "same length as pwc.vec" , ii goes from 0>mesh.ne-1

        element_dofs = fes.GetDofNrs(el) # tuple with the indices of the vector \psi2 corresponding to current element starting from zero 
        psi0_5 = np.zeros(6)

        # for k,dof in enumerate(element_dofs): #filling the vector psi0_5 with the required values from the \psi2 vector
        #     if k <=2:
        #         psi0_5[k]=levelset_p2.vec[dof]


        # # psi0_5[3]=levelset_p2.vec[element_dofs[5]]
        # # psi0_5[4]=levelset_p2.vec[element_dofs[4]]
        # # psi0_5[5]=levelset_p2.vec[element_dofs[3]]
        # psi0_5[3] =0.5*(psi0_5[0]+psi0_5[2])
        # psi0_5[4] =0.5*(psi0_5[1]+psi0_5[2])
        # psi0_5[5] =0.5*(psi0_5[0]+psi0_5[1])

           # print(dof)
        x_coords = np.zeros(6)
        y_coords = np.zeros(6)
        positive_count = np.sum(psi0_5 > EPS)
        negative_count = np.sum(psi0_5 < -EPS)
        zero_count = np.sum(abs(psi0_5) <EPS)
        #print("#+ , #-, #0",positive_count,negative_count,zero_count)
        print("Element Number:",ii)
        for j, v in enumerate(el.vertices):  #(el.vertices) is a tuple with three numbers corresppoinding to same indexing as in element_dof
            meshv = mesh[v]
            x_coords[j] = meshv.point[0]
            y_coords[j] = meshv.point[1]
           
        # Coordinates of the mid points
        x_coords[3] = 0.5 * (x_coords[0] + x_coords[2])
        x_coords[4] = 0.5 * (x_coords[1] + x_coords[2])
        x_coords[5] = 0.5 * (x_coords[0] + x_coords[1])
        
        y_coords[3] = 0.5 * (y_coords[0] + y_coords[2])
        y_coords[4] = 0.5 * (y_coords[1] + y_coords[2])
        y_coords[5] = 0.5 * (y_coords[0] + y_coords[1])

        mip0 = mesh(x_coords[0],y_coords[0])
        mip1 = mesh(x_coords[1],y_coords[1])
        mip2 = mesh(x_coords[2],y_coords[2])
        mip3 = mesh(x_coords[3],y_coords[3])
        mip4 = mesh(x_coords[4],y_coords[4])
        mip5 = mesh(x_coords[5],y_coords[5])

        PSI= CoefficientFunction(levelset_p2)
        psi0_5[0] = PSI(mip0)
        psi0_5[1] = PSI(mip1)
        psi0_5[2] = PSI(mip2)
        psi0_5[3] = PSI(mip3)
        psi0_5[4] = PSI(mip4)
        psi0_5[5] = PSI(mip5)

        
        s = GetVolumeFraction2(psi0_5, x_coords, y_coords, EPS)  # s... volume fraction of negative part of triangle
        print(",s=",s)
        nCutElems += 1
        func_p0.vec[ii] = val2 + s * (val1 - val2)  # == s*val1 + (1-s)*val2
