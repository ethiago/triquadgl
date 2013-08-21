#ifndef FASTMARCHING_H
#define FASTMARCHING_H

#include <QObject>
#include <QImage>
#include <queue>

class FastMarching : public QObject
{
    Q_OBJECT

    QImage m_src;

    bool **map;
    float **values;

public:
    explicit FastMarching(const QImage& src, QObject *parent = 0);

    void setSource(const QImage& src);

    QImage run();
    
private:

    void flood();
    bool isValid(int i, int j);
    void allocMats();
    void freeMats();
    float maxValues();
    

    class TexelStruct{
    public:
        int i,j;
        float value;
        TexelStruct(int _i, int _j, float v):i(_i),j(_j), value(v){}

        TexelStruct& operator=(const TexelStruct& t){
            i = t.i; j = t.j; value = t.value;
        }

        bool operator<(const TexelStruct& t)const{
            return value > t.value;
        }
    };

    std::priority_queue<TexelStruct> fila;
    
};

#endif // FASTMARCHING_H
