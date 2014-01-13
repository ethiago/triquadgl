#ifndef FASTMARCHING_H
#define FASTMARCHING_H

#include <QObject>
#include <QImage>
#include <queue>

class FastMarching : public QObject
{
    Q_OBJECT

    QImage m_src;

    bool **m_map;
    float **m_values;

    bool alreadyRun;

public:
    explicit FastMarching(const QImage& src, QObject *parent = 0);
    virtual ~FastMarching();

    const float ** values() const;

    float distanceTo(const FastMarching& fm);
    const QImage &getImage()const;
    QImage getMapImage()const;

    void run();

    void normalize(float t = -1);

    float maxValues()const;
    
private:

    void flood();
    bool isValid(int i, int j)const;
    void allocMap();
    void freeMap();

    

    class TexelStruct{
    public:
        int i,j;
        float value;
        TexelStruct(int _i, int _j, float v):i(_i),j(_j), value(v){}

        TexelStruct& operator=(const TexelStruct& t){
            i = t.i; j = t.j; value = t.value;
            return *this;
        }

        bool operator<(const TexelStruct& t)const{
            return value > t.value;
        }
    };

    std::priority_queue<TexelStruct> fila;
    
};

#endif // FASTMARCHING_H
